import os
from pathlib import Path
from typing import Dict, Optional
import requests

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

    def _resolve_rdrpcatch_path(self):
        """Automatically detect correct database path structure"""
        # Case 1: Direct rdrpcatch_dbs path
        if self.db_dir.name == "rdrpcatch_dbs":
            return self.db_dir

        # Case 2: Parent directory containing rdrpcatch_dbs
        candidate = self.db_dir / "rdrpcatch_dbs"
        if candidate.exists() and candidate.is_dir():
            return candidate

        # Case 3: Already contains hmm_dbs subdirectory
        hmm_check = self.db_dir / "hmm_dbs"
        if hmm_check.exists():
            return self.db_dir

        # Case 4: Fallback to original path
        return self.db_dir

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
        db_dir = self._resolve_rdrpcatch_path()
        for root,dirs,files in os.walk(db_dir):
            for name in dirs:
                if name == db_name:
                    for file in os.listdir(os.path.join(root, name)):
                        if file.endswith(".h3m"):
                            db_fn = file.rsplit(".", 1)[0]
                            db_path = os.path.join(root,name, db_fn)
                        else:
                            continue

        if not db_path:
            raise FileNotFoundError(f"{db_name} not found in {db_dir}")
        else:
            return db_path


    def fetch_mmseqs_db_path(self, db_name):
        """
        Fetches MMseqs database from the RdRpCATCH repository
        """
        if not os.path.exists(self.db_dir):
            raise FileNotFoundError(f"db_dir does not exist {self.db_dir}")

        db_dir = self._resolve_rdrpcatch_path()

        db_path = None
        for root,dirs,files in os.walk(db_dir):
            for dir in dirs:
                if dir == "mmseqs_dbs":
                    for dir_ in os.listdir(os.path.join(root, dir)):
                        if dir_ == db_name:
                            for file in os.listdir(os.path.join(root, dir, dir_)):
                                if file.endswith(".lookup"):
                                    db_fn = file.rsplit(".", 1)[0]
                                    db_path = os.path.join(root, dir, dir_, db_fn)
                                else:
                                    continue

        if not db_path:
            raise FileNotFoundError(f"{db_name} not found in {db_dir}")
        else:
            return db_path






class ZenodoDownloader:
    """Handles Zenodo database downloads using record IDs for version tracking"""

    def __init__(self, concept_doi: str, db_dir: Path):
        self.concept_doi = concept_doi
        self.db_dir = db_dir
        self.temp_dir = db_dir / "temp"
        self.lock_file = db_dir / ".lock"
        self.version_file = db_dir / "version.json"
        self._api_base = "https://zenodo.org/api/records"

        self.db_dir.mkdir(parents=True, exist_ok=True)

    def _get_record_id(self) -> str:
        """Extract numeric concept ID from Concept DOI"""
        return self.concept_doi.split(".")[-1]

    def _fetch_latest_metadata(self) -> Dict:
        """Retrieve latest version metadata from Zenodo API"""

        response = requests.get(f"{self._api_base}/{self._get_record_id()}")
        response.raise_for_status()
        return response.json()

    def get_latest_version_info(self) -> Dict:
        """Get version information from Zenodo metadata"""
        metadata = self._fetch_latest_metadata()
        return {
            "record_id": str(metadata["id"]),
            "doi": metadata["metadata"]["doi"],
            "created": metadata["metadata"]["publication_date"],
            "conceptdoi": self.concept_doi
        }

    def _get_tarball_info(self) -> Dict:
        """Find database tarball in Zenodo files"""

        metadata = self._fetch_latest_metadata()
        for file_info in metadata.get("files", []):
            if file_info["key"].endswith(".tar"):
                return {
                    "url": file_info["links"]["self"],
                    "checksum": file_info["checksum"],
                    "size": file_info["size"]
                }
        raise ValueError("No database tarball found in Zenodo record")


    def _verify_checksum(self, file_path: Path, expected: str) -> bool:
        """Validate file checksum (supports MD5 and SHA-256)"""
        import hashlib

        algorithm, _, expected_hash = expected.partition(":")
        hasher = hashlib.new(algorithm)

        with open(file_path, "rb") as f:
            while chunk := f.read(8192):
                hasher.update(chunk)

        return hasher.hexdigest() == expected_hash

    def extract_and_verify(self, tar_path: Path) -> None:
        """Safely extract tarball with proper directory structure handling"""
        import tarfile
        import tempfile
        import shutil


        with tempfile.TemporaryDirectory(dir=self.db_dir, prefix=".tmp_") as tmp_extract:
            tmp_extract = Path(tmp_extract)

            # Extract to temporary subdirectory
            with tarfile.open(tar_path, "r") as tar:
                # Validate all members first
                for member in tar.getmembers():
                    if ".." in member.path:
                        raise ValueError(f"Invalid path in archive: {member.path}")
                tar.extractall(tmp_extract)

            # Handle nested directory structure
            extracted_items = list(tmp_extract.iterdir())
            if len(extracted_items) == 1 and extracted_items[0].is_dir() and extracted_items[0].name == "rdrpcatch_dbs":
                # If archive contains single rdrpcatch_dbs directory, move contents up
                nested_dir = tmp_extract / "rdrpcatch_dbs"
                for item in nested_dir.iterdir():
                    shutil.move(str(item), str(tmp_extract))
                nested_dir.rmdir()


            # Prepare paths for atomic replacement
            target = self.db_dir / "rdrpcatch_dbs"
            backup = self.db_dir / "rdrpcatch_dbs.bak"

            # Atomic replacement sequence
            try:
                if backup.exists():
                    shutil.rmtree(backup)

                if target.exists():
                    target.rename(backup)

                tmp_extract.rename(target)

            finally:
                if backup.exists() and target.exists():
                    shutil.rmtree(backup)

    def needs_update(self) -> bool:
        """Check if local databases are outdated using record ID"""
        import json

        if not self.version_file.exists():
            return True

        try:
            with open(self.version_file, "r") as f:
                local_version = json.load(f)
            remote_version = self.get_latest_version_info()
            return remote_version["record_id"] != local_version["record_id"]
        except (json.JSONDecodeError, KeyError):
            return True


    def atomic_write_version(self, version_info: Dict) -> None:
        """Safely update version file with download timestamp"""
        import json
        import datetime

        temp_version = self.version_file.with_suffix(".tmp")

        # Add timestamp BEFORE writing to file
        version_info["downloaded"] = datetime.datetime.utcnow().isoformat()
         # Using ISO 8601 with timezone

        with open(temp_version, "w") as f:
            json.dump(version_info, f, indent=2)

        os.replace(temp_version, self.version_file)

    def get_current_version(self) -> Optional[Dict]:
        """Read installed database version info"""
        import json

        if self.version_file.exists():
            with open(self.version_file, "r") as f:
                return json.load(f)
        return None
