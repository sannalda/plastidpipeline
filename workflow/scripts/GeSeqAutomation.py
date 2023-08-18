import time
import os

from selenium.webdriver.chrome.service import Service as ChromeService
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import ElementClickInterceptedException
from webdriver_manager.chrome import ChromeDriverManager
#import chromedriver_binary 



##### Input Parameters
annotation_file_input = os.path.join(snakemake.config["workdir"],snakemake.input[0])
annotation_file_output = os.path.join(snakemake.config["workdir"],snakemake.output[0])


##### Config Parameters
GenomeShape = snakemake.config["GenomeShape"] 
SequenceSource = snakemake.config["SequenceSource"] 
AnnotateIR = snakemake.config["AnnotateIR"] 
AnnotateRPS12 = snakemake.config["AnnotateRPS12"] 
AnnotationChloe = snakemake.config["AnnotationChloe"] 
AnnotationMFannot = snakemake.config["AnnotationMFannot"] 
AnnotationRevision = snakemake.config["AnnotationRevision"] 
MPIMP_RefSet = snakemake.config["MPIMP_RefSet"] 
MultiGenBank = snakemake.config["MultiGenBank"] 



##### Selenium Webdriver
op = webdriver.ChromeOptions()
op.add_argument("--headless")
op.add_argument("--no-sandbox")
op.add_argument("--disable-dev-shm-usage")
driver = webdriver.Chrome(options=op)
driver.get("https://chlorobox.mpimp-golm.mpg.de/geseq.html")
time.sleep(5)



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
upload_fasta_input = upload_fasta_block.find_element(By.CLASS_NAME,"fl_head")
upload_fasta_input_file = upload_fasta_input.find_element(By.XPATH, "//input[@type='file']")
upload_fasta_input_file.send_keys(annotation_file_input)
#upload_fasta_input_file.send_keys("/Users/SJAnnaldasula/Documents/BGBM/Plastid/Am09/ChloroplastsCircular09Sequence1_AM0909.fasta")
#upload_fasta_input_file.send_keys(annotation_filename_standardized)
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

        
        
##### Right Column
right_column_elements = right_column.find_elements(By.CLASS_NAME, 'gs_panel')
output_options_block = right_column_elements[0]
actions_block = right_column_elements[1]
results_block = right_column_elements[2]

### Output options
assert(type(MultiGenBank) == bool)
       
if (MultiGenBank):
    if (not output_options_block.find_element(By.ID,"multigenbank_enabled").is_selected()):
        try:
            multigenbank_button = WebDriverWait(driver, 20).until(EC.element_to_be_clickable(output_options_block.find_element(By.ID,"multigenbank_enabled")))
            multigenbank_button.click()
        except ElementClickInterceptedException:
            print("Trying to click on the button again")
            driver.execute_script("arguments[0].click()", multigenbank_button)


### Actions
# Disclaimer
if (not actions_block.find_element(By.ID,"cb_disclaimer").is_selected()):
    actions_block.find_element(By.ID,"cb_disclaimer").click()
    
# Submit job
actions_block.find_element(By.CLASS_NAME,"gs_submit").click()
time.sleep(1)
submit_job_popup = driver.find_element(By.ID,"io_dialog")
job_title = submit_job_popup.find_element(By.CLASS_NAME,"input_string").get_attribute("value")
submit_job_popup.find_element(By.CLASS_NAME,"cms_button_ok").click()


# Waiting for job to run
time.sleep(1)
assert(job_title == results_block.find_element(By.CLASS_NAME,"gs_jobtitle").text.strip())

ignored_exceptions=(NoSuchElementException,StaleElementReferenceException,)

job_status = results_block.find_element(By.CLASS_NAME,"gs_jobstatus").text.strip()
while(job_status != 'Status: finished'):
    WebDriverWait(driver,5,ignored_exceptions=ignored_exceptions).until(EC.presence_of_element_located((By.CLASS_NAME,"gs_jobstatus")))
    time.sleep(1)
    job_status = results_block.find_element(By.CLASS_NAME,"gs_jobstatus").text.strip()
    print(job_status)
    
    
    
##### Downloading GenBank file
annotation_filename = results_block.find_element(By.XPATH,'//a[@data-gs-format="GenBank"]').get_attribute("data-gs_filename")
results_block.find_element(By.XPATH,'//a[@data-gs-format="GenBank"]').click()
time.sleep(2)
download_file_popup = driver.find_element(By.ID,"io_dialog")
download_file_popup.find_element(By.CLASS_NAME,"cms_button_download").click()
time.sleep(4)
driver.quit()        


##### Rename file
os.rename(annotation_filename,annotation_file_output)