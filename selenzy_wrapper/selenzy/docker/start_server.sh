#!/usr/bin/env bash
docker build -t sbc/selenzybase -f Dockerfile.base .
docker build --no-cache=true -t sbc/selenzy .
docker run -d -p 5555:5000 -e LD_LIBRARY_PATH='/opt/conda/bin/../lib' sbc/selenzy
# docker run -p 5555:5000 -e LD_LIBRARY_PATH='/opt/conda/bin/../lib' -it --entrypoint /bin/bash sbc/selenzy
# python /selenzyPro/flaskform.py -uploaddir /selenzyPro/uploads -datadir /selenzyPro/data -logdir /selenzyPro/log
