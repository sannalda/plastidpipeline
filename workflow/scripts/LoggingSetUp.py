import logging, os, time

sample = snakemake.params[0]

logger_sample = logging.getLogger("PlastidPipeline_"%snakemake.wildcards["sample"])
logger_sample.setLevel(logging.DEBUG)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
console_handler.setFormatter(formatter)
logger_sample.addHandler(console_handler)

time_start = datetime.now()
#logging_file = f"PlastidPipeline_{sample}_{time_start:%Y%m%d_%H%M%S}.log"
#log_file_path = copy.deepcopy(os.path.join(config["workdir"],"logs",logging_file))
file_handler = logging.FileHandler(snakemake.output[0])
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger_sample.addHandler(file_handler)


logger_sample.debug(snakemake.config)