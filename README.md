# mqCrunchAdmin


## to build the image run 
```
sudo docker build --no-cache . -t mqcruncher
sudo docker run -d -it --restart unless-stopped -v  $PWD/example/:/home/shiny/my_projects -v $PWD/example/table/:/home/shiny/table -v $PWD/example/annotations/:/home/shiny/annotations  -p 3839:3839 mqcruncher
```




