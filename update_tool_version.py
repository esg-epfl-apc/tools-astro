import re

def increment_version(file_path):
    version_pattern = r'(\d+\.\d+\.)(\d+)'

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        first_line = lines[0]
        match = re.search(version_pattern, first_line)

        if match:
            base_version = match.group(1)
            last_number = int(match.group(2))
            new_last_number = last_number + 1
            new_version = f"{base_version}{new_last_number}"
            updated_line = first_line.replace(match.group(0), new_version)
            lines[0] = updated_line

            with open(file_path, 'w') as file:
                file.writelines(lines)

    except FileNotFoundError:
        print(f"No archives file at: {file_path}")
    except Exception as e:
        print(f"Error while updating version: {e}")


if __name__ == "__main__":
    file_path = 'tools/archives/pyvo_integration/astronomical_archives.xml'

    increment_version(file_path)
