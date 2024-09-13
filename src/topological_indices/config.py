from pathlib import Path
import tomli

data_dir = (Path(__file__).parent / 'data').resolve()

with open((Path(__file__).parents[2]).resolve() / 'config.toml', "rb") as f:
    config = tomli.load(f)

nauty_path = Path(config['external']['nauty_path'])
viennarna_path = Path(config['external']['viennarna_path'])
