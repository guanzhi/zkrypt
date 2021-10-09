
all:
	gcc -Wall bn254.c bn254_params.c bn254_test.c -o bn254_test
	gcc -Wall bn254.c bn254_params.c plonk.c sha256.c plonk_test.c -o plonk_test

clean:
	rm -fr a.out
	rm -fr bn254_test
	rm -fr plonk_test

