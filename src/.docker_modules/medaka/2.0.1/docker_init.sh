    #!/bin/sh
# docker pull xgrand/medaka:2.0.1
docker build src/.docker_modules/medaka/2.0.1 -t 'xgrand/medaka:2.0.1'
docker push xgrand/medaka:2.0.1
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/medaka:2.0.1" --push src/.docker_modules/medaka/2.0.1
