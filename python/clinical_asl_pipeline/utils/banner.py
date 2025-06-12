# banner.py
import logging
def log_pipeline_banner(tool_version, logging=logging):
    logging.info("==============================================")
    logging.info(" ClinicalASL - Clinical Arterial Spin Labeling ")
    logging.info(f" Release: v{tool_version}                      ")
    logging.info(" Repository: https://github.com/JSIERO/ClinicalASL")
    logging.info(" Author: Jeroen Siero (UMC Utrecht, NL)        ")
    logging.info(" Contact: j.c.w.siero@umcutrecht.nl            ")
    logging.info("==============================================")