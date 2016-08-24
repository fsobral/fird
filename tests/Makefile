hs:
	$(MAKE) -C hs lib install

probability:
	$(MAKE) -C probability lib install

clean:
	rm -vf *~
	$(foreach i, $(shell ls -d */), $(MAKE) -C $(i) clean;)

.PHONY: hs probability clean
