from sqlalchemy import create_engine
import hashlib
import yaml


def get_sql_engine(config_name, config_file):
    with open(config_file,'r') as file:
        configs = yaml.safe_load(file)
    cfg = configs[config_name]
    user = cfg['user']
    password = cfg['password']
    host = cfg['host']
    port = cfg['port']
    database = cfg['database']

    engine = create_engine(f'mysql+mysqlconnector://{user}:{password}@{host}:{port}/{database}', echo=False)

    return engine


def get_seq_hash(sequence):
    return hashlib.md5(sequence.encode(encoding='utf-8')).hexdigest()


def collect_value(x, val):
    cur_dict = {k.split("=")[0]:k.split("=")[1] for k in x.split(";")}
    
    if val in cur_dict.keys():
        protein = cur_dict[val]
    else:
        protein = "."
    return protein