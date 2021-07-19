sudo docker run --name selenzy -v selenzyPro:/selenzy -p 5555:5000 -e LD_LIBRARY_PATH='/opt/conda/bin/../lib' -it selenzy python /selenzy/flaskform.py /selenzy/uploads /selenzy/data
