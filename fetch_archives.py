import pyvo


def fetch_archives():
    archive_list = pyvo.registry.search(servicetype="tap")

    with open('tools/archives/pyvo_integration/tool-data/astronomical_archives_gen.loc', 'w') as f, open('tools/archives/pyvo_integration/tool-data/astronomical_archives_gen.loc.sample', 'w') as f1:

        table_description = '#Table used for listing astronomical archives TAP service urls\n'
        column_description = '#<id>\t<display_name>\t<value>\n'

        f.write(table_description)
        f.write(column_description)

        f1.write(table_description)
        f1.write(column_description)

        archive = archive_list[0]
        archive_name = archive.res_title
        access_url = archive.access_url
        
        f1.write(f'0\t{archive_name}\t{access_url}\n')
        
        for i, archive in enumerate(archive_list):
            try:
                archive_name = archive.res_title
                access_url = archive.access_url
                f.write(f'{i}\t{archive_name}\t{access_url}\n')
            except Exception as e:
                pass


if __name__ == '__main__':
    fetch_archives()

