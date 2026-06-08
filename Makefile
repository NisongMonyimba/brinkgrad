.PHONY: reproduce test clean

reproduce:
	docker run --rm \
		-v $(PWD):/root/brinkgrad \
		-e OMP_NUM_THREADS=1 \
		-e MPLBACKEND=Agg \
		dolfinx/dolfinx:v0.7.3 \
		bash -c "cd /root/brinkgrad && pip install -e . --quiet && bash run_all.sh"

test:
	docker run --rm \
		-v $(PWD):/root/brinkgrad \
		-e OMP_NUM_THREADS=1 \
		dolfinx/dolfinx:v0.7.3 \
		bash -c "cd /root/brinkgrad && pip install -e . --quiet && python -m pytest tests/ -v"

clean:
	rm -rf figures/*.pdf figures/*.png figures/*.csv
