# mqCrunchAdmin


## to build the image run 
```
sudo docker build . -t mqctest
```

## to deploy the image as a container on Linux Computing3
```
sh deploy_as_container.sh
```
or

```
docker-compose up -d --build

```
## to test if the container is up and running
```
sudo docker ps
```

## to register the image in the dockerHub
```
docker login
docker tag mqctest:latest nava20ir/mqctest:latest
docker push nava20ir/mqctest:latest

```


