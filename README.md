# mqCrunchAdmin

How this repository was deconvoluted from the docker container can be foud at:
https://collab.dvb.bayern/spaces/TUMbaybioms/pages/1721605749/How+to+deconvolure+one+R+package+from+a+docker+i.e+MQCrunch

## Getting started
all the volumes needs to be mounted before deployment

https://collab.dvb.bayern/spaces/TUMbaybioms/pages/1486370996/Mounting+netapp+drives+on+linux_computing+3


## to build the image run 
```
sudo docker build . -t mqctest
```

## to deploy the image as a container on Linux Computing3
```
sh deploy_as_container.sh
```

## to test if the container is up and running
```
sudo docker ps
```


