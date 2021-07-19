# Docker file that installs docker container for Selenzy
#
# rename this file to "Dockerfile"
# build with: "sudo docker build -t selenzy ."

FROM sbc/selenzybase 

# To be replaced by a git clone
RUN wget http://130.88.113.226/selenzy/selenzy.tar.gz
RUN tar -xzvf selenzy.tar.gz

ENTRYPOINT ["python"]

CMD ["/selenzyPro/flaskform.py", "-uploaddir", "/selenzyPro/uploads", "-datadir", "/selenzyPro/data", "-logdir", "/selenzyPro/log" ]

EXPOSE 5000


