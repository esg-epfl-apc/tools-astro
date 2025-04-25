import os
import json

from minio import Minio

# from minio.error import S3Error


class CRbeamCache(object):
    credentials_env_var = "S3_CREDENTIALS"
    bucket_name = "crbeam"
    chunk_size = 32 * 1024

    # default credentials
    endpoint = "play.min.io"
    access_key = "Q3AM3UQ867SPQQA43P2F"
    secret_key = "zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG"

    def __init__(self) -> None:
        credentials = os.getenv(self.credentials_env_var)
        if credentials is not None:
            credentials = json.loads(credentials)
            self.endpoint = credentials["endpoint"]
            self.access_key = credentials["access_key"]
            self.secret_key = credentials["secret_key"]

        self.client = Minio(
            self.endpoint,
            access_key=self.access_key,
            secret_key=self.secret_key,
        )

        # Make bucket if not exist.
        if not self.client.bucket_exists(self.bucket_name):
            self.client.make_bucket(self.bucket_name)

    def save(self, obj_name, file_path, append_mode=True, **params):
        if append_mode:
            idx = 0
            for _ in self.client.list_objects(
                self.bucket_name, prefix=obj_name
            ):
                idx += 1
            obj_path = obj_name + f"-{idx}"
        else:
            obj_path = obj_name
        self.client.fput_object(
            self.bucket_name, obj_path, file_path, metadata=params
        )
        return obj_path

    def get_cache_size(self, prefix):
        size = 0
        for obj in self.client.list_objects(self.bucket_name, prefix=prefix):
            size += int(
                obj.object_name[len(prefix):].split("-")[0]
            )  # todo: load number of particles from metadata
        return size

    def load_file(self, output_path, obj_name):
        try:
            response = self.client.get_object(self.bucket_name, obj_name)
        except Exception:
            raise ValueError("object not found")
        try:
            # Read data from response.
            with open(output_path, "wb") as file:
                for d in response.stream(self.chunk_size):
                    file.write(d)
        finally:
            response.close()
            response.release_conn()

    def load_results(self, output_path, prefix, skip_paths=[]):
        for obj in self.client.list_objects(self.bucket_name, prefix=prefix):
            print("found", obj.object_name)
            if obj.object_name not in skip_paths:
                # Get data of an object.
                response = self.client.get_object(
                    self.bucket_name, obj.object_name
                )
                try:
                    # Append the object data to a local file
                    with open(output_path, "ab") as file:
                        # todo: cut first line
                        for d in response.stream(self.chunk_size):
                            file.write(d)
                    # Read data from response.
                finally:
                    response.close()
                    response.release_conn()

    def detete_results(self, prefix):
        for obj in self.client.list_objects(self.bucket_name, prefix=prefix):
            self.client.remove_object(self.bucket_name, obj.object_name)


if __name__ == "__main__":
    # test code
    obj_name_prefix = "photon_z0.1_E_1e+07_1e+13_B1e-15_L0.05-5_N"
    n_particles = 10000
    root_path = "/Users/ok/git/mcray/bin/"
    particle = "photon"
    obj_name = obj_name_prefix + str(n_particles)
    file_dir = root_path + obj_name + "/z0"
    file_path = file_dir + "/" + particle
    c = CRbeamCache()
    c.save(obj_name, file_path, n_particles=n_particles)
    print(file_path + " sent to s3")
    loaded_file_path = file_path + ".loaded"
    c.load_results(loaded_file_path, obj_name)
    print(loaded_file_path + " loaded from s3")
