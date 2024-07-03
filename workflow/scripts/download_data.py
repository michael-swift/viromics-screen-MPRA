import json
import requests
import os
import sys

def download_file(url, output_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
    else:
        print(f"Failed to download {url}")

def parse_and_download(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    for sample in data['samples']:
        url = sample['url']
        output_path = os.path.join("data", f"{sample['name']}.fq.gz")
        if not os.path.exists(output_path):
            print(f"Downloading {url} to {output_path}")
            download_file(url, output_path)
        else:
            print(f"File {output_path} already exists, skipping download.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python download_data.py <json_file>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    parse_and_download(json_file)
