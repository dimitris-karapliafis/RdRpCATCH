import os

class db_fetcher:

    def __init__(self,db_dir):
        self.db_dir = db_dir
        self.version_file = os.path.join(db_dir, "version.json")
        self.custom_db_dir = os.path.join(db_dir, "custom_dbs")

    def _get_db_version(self):
        """Get database version information"""
        import json
        
        if os.path.exists(self.version_file):
            with open(self.version_file) as f:
                return json.loads(f.read())
        return {}

    def _save_db_version(self, version_info):
        """Save database version information"""
        import json
        
        os.makedirs(os.path.dirname(self.version_file), exist_ok=True)
        with open(self.version_file, 'w') as f:
            json.dump(version_info, f, indent=2)

    def check_db_updates(self):
        """Check if database updates are available"""
        current_version = self._get_db_version()
        # TODO: Implement version checking against remote repository
        # For now just return the current version
        return current_version

    def add_custom_db(self, db_path, db_name=None):
        """Add a custom database (MSA or pHMM file) to the custom_dbs directory"""
        import shutil
        import datetime
        
        if not os.path.exists(self.custom_db_dir):
            os.makedirs(self.custom_db_dir)

        if db_name is None:
            db_name = os.path.basename(db_path)

        target_path = os.path.join(self.custom_db_dir, db_name)
        
        # Copy the database file
        if os.path.isfile(db_path):
            shutil.copy2(db_path, target_path)
        elif os.path.isdir(db_path):
            if os.path.exists(target_path):
                shutil.rmtree(target_path)
            shutil.copytree(db_path, target_path)

        # Update version info
        version_info = self._get_db_version()
        version_info.setdefault('custom_dbs', {})
        version_info['custom_dbs'][db_name] = {
            'added': datetime.datetime.now().isoformat(),
            'path': target_path
        }
        self._save_db_version(version_info)

    def fetch_hmm_db_path(self, db_name):
        """
        Fetches HMM database from the RdRpCATCH repository or custom databases
        """
        if not os.path.exists(self.db_dir):
            raise FileNotFoundError(f"db_dir does not exist {self.db_dir}")

        # First check custom databases
        if os.path.exists(self.custom_db_dir):
            custom_path = os.path.join(self.custom_db_dir, db_name)
            if os.path.exists(custom_path):
                if os.path.isfile(custom_path) and custom_path.endswith(('.h3m', '.hmm')):
                    return os.path.splitext(custom_path)[0]
                elif os.path.isdir(custom_path):
                    for file in os.listdir(custom_path):
                        if file.endswith(('.h3m', '.hmm')):
                            return os.path.splitext(os.path.join(custom_path, file))[0]

        # Then check standard databases
        db_path = None
        for dir in os.listdir(self.db_dir):
            if dir == "hmm_dbs":
                for dir_ in os.listdir(os.path.join(self.db_dir, dir)):
                    if dir_ == db_name:
                        for file in os.listdir(os.path.join(self.db_dir, dir, dir_)):
                            if file.endswith(".h3m"):
                                db_fn = file.rsplit(".",1)[0]
                                db_path = os.path.join(self.db_dir, dir, dir_, db_fn)
                            else:
                                continue

        if not db_path:
            raise FileNotFoundError(f"{db_name} not found in {self.db_dir}")
        else:
            return db_path


    def fetch_mmseqs_db_path(self, db_name):
        """
        Fetches MMseqs database from the RdRpCATCH repository
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




