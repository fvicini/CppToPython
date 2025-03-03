.PHONY: create distclean clean run run-bash

## Clean
clean:
	@docker container prune
	@docker rmi $$(docker images --filter "dangling=true" -q --no-trunc)

## Create
create:
	@docker build -f Dockerfile --target rbm-base . -t rbm-base:1.0.0

## Run
run:
	@docker run --rm -it -p 8080:8080 rbm-base:1.0.0

run-bash:
	@docker run --rm -it -p 8080:8080 rbm-base:1.0.0 bash

## Distclean
distclean:
	@docker container prune
	@docker rmi rbm-base:1.0.0
	@docker rmi $$(docker images --filter "dangling=true" -q --no-trunc)
	@docker buildx prune -f
