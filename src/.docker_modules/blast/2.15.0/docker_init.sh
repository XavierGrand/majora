#!/bin/sh
# docker pull xgrand/blast:2.15.0
# docker build src/.docker_modules/blast/2.15.0 -t 'xgrand/blast:2.15.0'
# docker push xgrand/blast:2.15.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/blast:2.15.0" --push src/.docker_modules/blast/2.15

# docker pull xgrand/blast/2.15
docker build src/.docker_modules/blast/2.15 -t 'xgrand/blast:2.15.0'
docker push xgrand/blast:2.15.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/blast:2.15.0" --push src/.docker_modules/blast/2.15