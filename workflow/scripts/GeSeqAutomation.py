import time, os, logging

from selenium.webdriver.chrome.service import Service as ChromeService
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import ElementClickInterceptedException
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service as ChromiumService
from webdriver_manager.chrome import ChromeDriverManager
from webdriver_manager.core.os_manager import ChromeType
from selenium.webdriver import ActionChains

def setup_logger(log_file_path):
    loggerPP = logging.getLogger("PlastidPipeline_%s" %snakemake.wildcards["sample"])
    loggerPP.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    console_handler.setFormatter(formatter)
    loggerPP.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    loggerPP.addHandler(file_handler)

    return loggerPP

class GeSeqError(Exception):
    pass

class DownloadingError(Exception):
    pass

##### Input Parameters
loggerPP = setup_logger(os.path.join(snakemake.config["workdir"],snakemake.log[0]))
annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[0])
annotation_file_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])




##### Config Parameters
GenomeShape = snakemake.config["GeSeq"]["GenomeShape"] 
SequenceSource = snakemake.config["GeSeq"]["SequenceSource"] 
AnnotateIR = snakemake.config["GeSeq"]["AnnotateIR"] 
AnnotateRPS12 = snakemake.config["GeSeq"]["AnnotateRPS12"] 
AnnotationChloe = snakemake.config["GeSeq"]["AnnotationChloe"] 
ChloeAnnotateCDS = snakemake.config["GeSeq"]["ChloeAnnotateCDS"]
ChloeAnnotateTRNA = snakemake.config["GeSeq"]["ChloeAnnotateTRNA"]  
ChloeAnnotateRRNA = snakemake.config["GeSeq"]["ChloeAnnotateRRNA"] 
AnnotationMFannot = snakemake.config["GeSeq"]["AnnotationMFannot"] 
AnnotationRevision = snakemake.config["GeSeq"]["AnnotationRevision"] 
MPIMP_RefSet = snakemake.config["GeSeq"]["MPIMP_RefSet"] 
MultiGenBank = snakemake.config["GeSeq"]["MultiGenBank"] 



##### Selenium Webdriver
#options = Options()
#options.BinaryLocation = "/usr/bin/chromium-browser" 

op = webdriver.ChromeOptions()
op.add_argument("--headless")
op.add_argument('--ignore-certificate-errors')
op.add_argument("--no-sandbox")
op.add_argument("--disable-dev-shm-usage")
op.add_argument("--start-maximized")
op.add_argument("--disable-gpu");

### Download location
prefs = {"download.default_directory": snakemake.config["workdir"]};
op.add_experimental_option("prefs", prefs);

############# Code from "https://github.com/SergeyPirogov/webdriver_manager/issues/664#issuecomment-2247178221" to fix "OSError: [Errno 8] Exec format error:" die to THIRD_PARTY_NOTICES.chromedriver being added in July 2024 
driver_path = ChromeDriverManager(chrome_type=ChromeType.CHROMIUM).install()
if driver_path:
    driver_name = driver_path.split('/')[-1]
    if driver_name!="chromedriver":
        driver_path = "/".join(driver_path.split('/')[:-1]+["chromedriver"])
        os.chmod(driver_path, 0o755)
#######################################

Service = ChromiumService(driver_path)
driver = webdriver.Chrome(options=op,service=Service)
#driver = webdriver.Chrome(options=op)
driver.get("https://chlorobox.mpimp-golm.mpg.de/geseq.html")
time.sleep(15)

##### Initialize tool

elements = driver.find_elements(By.CLASS_NAME, 'x4_column')
left_column = elements[0]
middle_column = elements[1]
right_column = elements[2]



##### Left column 

left_column_elements = left_column.find_elements(By.CLASS_NAME, 'gs_panel')
upload_fasta_block = left_column_elements[0]
annotation_options = left_column_elements[1]

### Fasta files to annotate
# Uploading fasta files to annotate
upload_fasta_input = upload_fasta_block.find_element(By.CLASS_NAME,"fl_head")####
upload_fasta_input_file = upload_fasta_input.find_element(By.XPATH, "//input[@type='file']")
upload_fasta_input_file.send_keys(annotation_file_input)
time.sleep(5)

# Choose circular or linear genome
assert(GenomeShape in ["Linear","Circular"])

if (GenomeShape == "Linear"):
    upload_fasta_block.find_element(By.ID,"ogd_shape_linear").click()
else:
    upload_fasta_block.find_element(By.ID,"ogd_shape_circular").click()
    
# Sequence source
assert(SequenceSource in ["Land","Algae","Mito"])

if (SequenceSource == "Land"):
    upload_fasta_block.find_element(By.ID,"gssrctype_chloro").click()
elif (SequenceSource == "Algae"):
    upload_fasta_block.find_element(By.ID,"gssrctype_chloro_algae").click()
else:
    upload_fasta_block.find_element(By.ID,"gssrctype_mito").click()
# Note: need to look at functionality to algae and mitochondria? (like some boxes are unable to be ticked)

# Annotation Options
assert(type(AnnotateIR) == bool and type(AnnotateRPS12) == bool)

if (AnnotateIR):
    if (not upload_fasta_block.find_element(By.ID,"annotate_repeats").is_selected()):
        upload_fasta_block.find_element(By.ID,"annotate_repeats").click()

if (AnnotateRPS12):
    if (not upload_fasta_block.find_element(By.ID,"rps12splicing").is_selected()):
        upload_fasta_block.find_element(By.ID,"rps12splicing").click()


# Annotation Support
assert(type(AnnotationChloe) == bool and type(AnnotationMFannot) == bool)

if (AnnotationChloe):
    if (not upload_fasta_block.find_element(By.ID,"support_anno_chloe").is_selected()):
        upload_fasta_block.find_element(By.ID,"support_anno_chloe").click()
#if (AnnotationMFannot):
#    if (not upload_fasta_block.find_element(By.ID,"support_anno_mfannot").is_selected()):
#        upload_fasta_block.find_element(By.ID,"support_anno_mfannot").click()
# Note: MFannot not really working (on GeSeq side) - or more of I cannot seem to find to parameters to enable it

# Annotation revision
assert(AnnotationRevision in ["Best","All"])

if (AnnotationRevision == "Best"):
    upload_fasta_block.find_element(By.ID,"dropmode_keepbest").click()
else:
    upload_fasta_block.find_element(By.ID,"dropmode_keepall").click()
    
    
    
##### Middle Column 
middle_column_elements = middle_column.find_elements(By.CLASS_NAME, 'gs_panel')
blat_ref_seqs_block = middle_column_elements[0]
third_party_annotation_block = middle_column_elements[1]

### BLAT Reference Sequences
# MPI-MP Reference Set
assert(type(MPIMP_RefSet) == bool)
       
if (MPIMP_RefSet):
    if (not blat_ref_seqs_block.find_element(By.ID,"mpimpchlororefsetenabled").is_selected()):
        blat_ref_seqs_block.find_element(By.ID,"mpimpchlororefsetenabled").click()

### 3rd Party Stand-Alone Annotators
# Chloe
assert(type(ChloeAnnotateCDS) == bool)
assert(type(ChloeAnnotateTRNA) == bool)
assert(type(ChloeAnnotateRRNA) == bool)

if (AnnotationChloe):
    if (ChloeAnnotateCDS):
        if (not third_party_annotation_block.find_element(By.ID,"chloe_annotate_cds").is_selected()):
            third_party_annotation_block.find_element(By.ID,"chloe_annotate_cds").click()

    if (ChloeAnnotateTRNA):
        if (not third_party_annotation_block.find_element(By.ID,"chloe_annotate_trna").is_selected()):
            third_party_annotation_block.find_element(By.ID,"chloe_annotate_trna").click()

    if (ChloeAnnotateRRNA):
        if (not third_party_annotation_block.find_element(By.ID,"chloe_annotate_rrna").is_selected()):
            third_party_annotation_block.find_element(By.ID,"chloe_annotate_rrna").click()


##### Right Column
right_column_elements = right_column.find_elements(By.CLASS_NAME, 'gs_panel')
output_options_block = right_column_elements[0]
actions_block = right_column_elements[1]
results_block = right_column_elements[2]

### Output options
assert(type(MultiGenBank) == bool)
       
if (MultiGenBank):
    if (not output_options_block.find_element(By.ID,"multigenbank_enabled").is_selected()):
        multigenbank_element = output_options_block.find_element(By.ID,"multigenbank_enabled")
        actions = ActionChains(driver)
        actions.move_to_element(multigenbank_element).click().perform()

        #try:
        #    multigenbank_button = WebDriverWait(driver, 20).until(EC.element_to_be_clickable(output_options_block.find_element(By.ID,"multigenbank_enabled")))
        #    multigenbank_button.click()
        #except ElementClickInterceptedException:
        #    print("Trying to click on the button again")
        #    driver.execute_script("arguments[0].click()", multigenbank_button)


### Actions
# Disclaimer
if (not actions_block.find_element(By.ID,"cb_disclaimer").is_selected()):
    actions_block.find_element(By.ID,"cb_disclaimer").click()
    
# Submit job
actions_block.find_element(By.CLASS_NAME,"gs_submit").click()
time.sleep(5)
submit_job_popup = driver.find_element(By.ID,"io_dialog")
job_name_element = submit_job_popup.find_element(By.CLASS_NAME,"input_string")
#driver.execute_script("arguments[0].setAttribute('value',arguments[1])",job_name_element,snakemake.wildcards["sample"])
driver.execute_script("arguments[0].setAttribute('value',arguments[1])",job_name_element,"AutomatedScript")
loggerPP.debug("Submmited job to GeSeq for sample %s..." %snakemake.wildcards["sample"])
time.sleep(2)
job_title = submit_job_popup.find_element(By.CLASS_NAME,"input_string").get_attribute('value')
submit_job_popup.find_element(By.CLASS_NAME,"cms_button_ok").click()

# Waiting for job to run
time.sleep(1)
assert(job_title == results_block.find_element(By.CLASS_NAME,"gs_jobtitle").text.strip())

ignored_exceptions=(NoSuchElementException,StaleElementReferenceException,)

start_time = time.time()
job_status = results_block.find_element(By.CLASS_NAME,"gs_jobstatus").text.strip()
while(job_status != 'Status: finished'):
    #WebDriverWait(driver,5,ignored_exceptions=ignored_exceptions).until(EC.presence_of_element_located((By.CLASS_NAME,"gs_jobstatus")))
    time.sleep(15)
    try:
        job_status = results_block.find_element(By.CLASS_NAME,"gs_jobstatus").text.strip()
    except StaleElementReferenceException as e:
        pass
    loggerPP.debug("\t"+job_status)
    curr_time = time.time()
    try:
        if (curr_time - start_time > 1000):
            raise GeSeqError("GeSeq took too long, most likely because could not access the server/website. Please try again later by rerunning the script. Exiting...")
    except GeSeqError as e:
        loggerPP.error("GESeqError error detected: %s", e)
        break

time.sleep(5)
##### Downloading GenBank file
annotation_filename = results_block.find_element(By.XPATH,'//a[@data-gs-format="GenBank"]').get_attribute("data-gs_filename")
results_block.find_element(By.XPATH,'//a[@data-gs-format="GenBank"]').click()
time.sleep(10)
download_file_popup = driver.find_element(By.ID,"io_dialog")
#time.sleep(2)
download_file_popup.find_element(By.CLASS_NAME,"cms_button_download").click()
start_time = time.time()

while not os.path.exists(os.path.join(snakemake.config["workdir"],annotation_filename)):
    loggerPP.debug("\tAnnotation completed. Wating for download to complete...")
    time.sleep(10)
    curr_time = time.time()
    try:
        if (curr_time - start_time > 500):
            raise DownloadingError("Downloading took too long, please try again later by rerunning the script. Exiting...")
    except DownloadingError as e:
        loggerPP.error("DownloadingError error detected: %s", e)
        break
    
time.sleep(10)
driver.quit()  
if not os.path.exists(os.path.join(snakemake.config["workdir"],annotation_filename)):
    loggerPP.info("Downloading of annotated file failed. Location:", os.getcwd(),"- Filename:", annotation_filename)
else:
    loggerPP.info("Annotation and download complete.")      



##### Rename file
os.rename(annotation_filename,annotation_file_output)
