#!/bin/sh
docker pull xgrand/rlollipop:1.1
# docker build src/.docker_modules/r-majora -t 'xgrand/rlollipop:1.1'
# docker push xgrand/rlollipop:1.1
docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/rlollipop:1.1" --push src/.docker_modules/r-majora