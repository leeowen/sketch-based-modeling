{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 import requests\
import os\
import posixpath\
\
url_template=r\'92http://cdn.gea.esac.esa.int/Gaia/tgas_source/csv/TgasSource_000-000-\{index:03d\}.csv.gz\'92\
\
def main():\
	for index in xrange(16):\
		url=url_template.format(index=index)\
		filename=posixpath.basename(url)\
		if os.path.isfile(filename):\
			print \'93skipped\'94+filename\
			continue\
		print \'93Downloading\'94, filename\
		response=requests.get(url,stream=True)\
		with open(filename,\'92wb\'92) as fp:\
			for chunk in response.iter_content(1024):\
				fp.write(chunk)\
\
if __name__==\'91__main()__\'92:\
	main()}