APPID ?= wallen-Autodock-Vina-1.2.3

list:
	tapis apps search --id start wallen

deploy:
	tapis apps deploy -W . -F app.json --ini project.ini --no-build --no-push

update:
	tapis apps update -F app.json --ini project.ini ${APPID}


build:
	docker build -t wjallen/autodock_vina:1.2.3 .

run:
	docker run --rm -it wjallen/autodock_vina:1.2.3 /bin/bash

push:
	docker push wjallen/autodock_vina:1.2.3

