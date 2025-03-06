import os



class db_fetcher:

    def __init__(self,db_dir):
        self.db_dir = db_dir


    def fetch_hmm_db_path(self, db_name):
        """
        Fetches HMM database from the colabscanner repository
        """
        if not os.path.exists(self.db_dir):
            raise FileNotFoundError(f"db_dir does not exist {self.db_dir}")

        db_path = None
        for dir in os.listdir(self.db_dir):

            if dir == "hmm_dbs":
                for dir_ in os.listdir(os.path.join(self.db_dir, dir)):
                    if dir_ == db_name:
                        for file in os.listdir(os.path.join(self.db_dir, dir, dir_)):
                            if file.endswith(".h3m") :
                                db_fn = file.rsplit(".",1)[0]
                                db_path = os.path.join(self.db_dir, dir,dir_, db_fn)
                            else:
                                continue

        if not db_path:
            raise FileNotFoundError(f"{db_name} not found in {self.db_dir}")
        else:
            return db_path


    def fetch_mmseqs_db_path(self, db_name):
        """
        Fetches MMseqs database from the colabscanner repository
        """
        if not os.path.exists(self.db_dir):
            raise FileNotFoundError(f"db_dir does not exist {self.db_dir}")

        db_path = None
        for dir in os.listdir(self.db_dir):

            if dir == "mmseqs_dbs":
                for dir_ in os.listdir(os.path.join(self.db_dir, dir)):
                    if dir_ == db_name:
                        for file in os.listdir(os.path.join(self.db_dir, dir, dir_)):
                            if file.endswith(".lookup"):
                                db_fn = file.rsplit(".",1)[0]
                                db_path = os.path.join(self.db_dir, dir,dir_, db_fn)
                            else:
                                continue

        if not db_path:
            raise FileNotFoundError(f"{db_name} not found in {self.db_dir}")
        else:
            return db_path


class db_downloader:

    def __init__(self,destination_dir):
        self.base_url = "https://zenodo.org/records/14936013/files/colabscan_dbs.tar?download=1"
        self.destination_dir = destination_dir

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
        os.system(f"wget -O {self.destination_dir}/colabscan_dbs.tar {self.base_url}")


    def extract_db(self):
        '''
        Extracts the hmm databases
        :return:
        '''
        print(f"Extracting hmm databases to {self.destination_dir}")
        os.system(f"tar -xvf {self.destination_dir}/colabscan_dbs.tar -C {self.destination_dir}")
        print(f"Extraction complete")

    def del_tar(self):
        '''
        Deletes the tar file
        :return:
        '''
        print(f"Deleting tar file")
        os.system(f"rm {self.destination_dir}/colabscan_dbs.tar")
        print(f"Deletion complete")




