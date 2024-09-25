from pathlib import Path
import yaml

output_folder = Path(config["output_folder"])
loom_files = Path(config["loom_files"])

loom_file_list = []
for p in loom_files.iterdir():
    if p.is_file() and p.suffix == ".loom":
        loom_file_list.append(p)


output_files = [output_folder / 'pools' / 'pool_(0,0).csv', output_folder / 'pools' / 'pool_(1,0).csv', output_folder / 'pools' / 'pool_(2,0).csv', output_folder / 'pools' / 'pool_(0,1).csv', output_folder / 'pools' / 'pool_(1,1).csv', output_folder / 'pools' / 'pool_(3,1).csv']

for file in loom_file_list:
    output_files.append(output_folder / 'temp' / f'{file.stem}_split1.loom')
    output_files.append(output_folder / 'temp' / f'{file.stem}_split2.loom')
    
def get_variant_list(pool_ID):
    file = output_folder / f"pool_{pool_ID}.txt"
    with open(file, 'r') as f:
        return [line.strip() for line in f]