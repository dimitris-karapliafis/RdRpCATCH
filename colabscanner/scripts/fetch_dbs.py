import os

class db_fetcher:
    def __init__(self, db_dir, db_name):
        self.db_dir = db_dir
        self.db_name = db_name

    def fetch_db_path(self):
        """
        Fetches the database from the colabscanner github repository
        """
        if not os.path.exists(self.db_dir):
            print(f"db_dir does not exist {self.db_dir}")
        db_path = None
        print(os.listdir(self.db_dir))
        for dir in os.listdir(self.db_dir):

            if dir == self.db_name:
                print(f"{self.db_name} found in {self.db_dir}")
                for file in os.listdir(os.path.join(self.db_dir, dir)):
                    if file.endswith(".hmm") or file.endswith(".h3m" or file.endswith(".h3i")) or file.endswith(".h3f") or file.endswith(".h3p"):
                        db_fn = file.rsplit(".",1)[0]
                        db_path = os.path.join(self.db_dir, dir, db_fn)
                    else:
                        continue


        if not db_path:
            print(f"{self.db_name} not found in {self.db_dir}")
        else:
            return db_path


class db_downloader:

    def __init__(self,destination_dir):
        self.base_url = "https://zenodo.org/records/14358349/files/hmm_dbs.tar?download=1"
        self.destination_dir = destination_dir
        self.filenames = []

    def download_db(self):
        '''
        Downloads the hmm databases from Zenodo
        :return:
        '''

        if not os.path.exists(self.destination_dir):
            os.makedirs(self.destination_dir)
        else:
            print(f"{self.destination_dir} already exists")

        print(f"Downloading hmm databases to {self.destination_dir}")
        os.system(f"wget -O {self.destination_dir}/hmm_dbs.tar {self.base_url}")


    def extract_db(self):
        '''
        Extracts the hmm databases
        :return:
        '''
        print(f"Extracting hmm databases to {self.destination_dir}")
        os.system(f"tar -xvf {self.destination_dir}/hmm_dbs.tar -C {self.destination_dir}")
        print(f"Extraction complete")

    def del_tar(self):
        '''
        Deletes the tar file
        :return:
        '''
        print(f"Deleting tar file")
        os.system(f"rm {self.destination_dir}/hmm_dbs.tar")
        print(f"Deletion complete")




