#!/bin/sh
# docker pull xgrand/vcflib:1.0.14
# docker build src/.docker_modules/vcflib/1.0.14  -t 'xgrand/vcflib:1.0.14'
# docker push xgrand/vcflib:1.0.14
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/vcflib:1.0.14" --push src/.docker_modules/vcflib/1.0.14 

# docker pull xgrand/seqkit:2.8.2
docker build src/.docker_modules/vcflib/1.0.14 -t 'xgrand/vcflib:1.0.14'
docker push xgrand/vcflib:1.0.14
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/vcflib:1.0.14" --push src/.docker_modules/vcflib/1.0.14 