# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 10:23:09 2018

@author: Jack
"""
import os
import sys
import smtplib
from email.mime.text import MIMEText

sbj = sys.argv[1]
anly = sys.argv[2]
test_user = sys.argv[3]

if anly.endswith('/'):
    anly = anly[:-1]

project = os.path.basename(anly)

msg = {}

if sbj == "Yes":
    title = "Testing SUCCESS! " + project
    body = "Hi, Guys\n\nThe project " + project + " has run successfully, Please check EVERYTHING and deliver the data.\n\n Path:\n\n" + anly + "\n\nThanks\n\nCCRSF_IFX"
else:
    title = "Testing Failed! " + project
    body = "Hi, Guys\n\nThe project " + project + " is failed, Please check the log file and try again.\n\n Path:\n\n" + anly + "\n\nThanks\n\nCCRSF_IFX"    


msg = MIMEText(body)
msg['Subject'] = title
me = "CCRSF_IFX@nih.gov"
msg['From'] = me

#you = ['keyur.talsania@nih.gov','xiongfong.chen2@nih.gov', 'vicky.chen@nih.gov', 'tsai-wei.shen@nih.gov', 'zhaoyong@mail.nih.gov', 'paul.schaughency@nih.gov', 'sulbha.choudhari@nih.gov', 'brittany.dulek@nih.gov', 'julie.park@nih.gov']
#msg['To'] = "keyur.talsania@nih.gov, xiongfong.chen2@nih.gov, vicky.chen@nih.gov, tsai-wei.shen@nih.gov, zhaoyong@mail.nih.gov, paul.schaughency@nih.gov, sulbha.choudhari@nih.gov, brittany.dulek@nih.gov, julie.park@nih.gov"
you = [test_user]
msg['To'] = test_user

s = smtplib.SMTP('localhost')
s.sendmail(me, you, msg.as_string())
s.quit()
