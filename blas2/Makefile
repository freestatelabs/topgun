
slow:
	gcc -O3 -march=native blas2_slow.c -o slow.o 
	./slow.o

v1:
	gcc -O3 -march=native -funroll-loops \
	-Rpass-missed=loop-vectorize -Rpass=loop-vectorize \
	-Rpass-analysis=loop-vectorize \
	-framework Accelerate \
	blas2_v1.c -o v1.o 
	./v1.o

v1_float:
	gcc -O3 -march=native -funroll-loops \
	-Rpass-missed=loop-vectorize -Rpass=loop-vectorize \
	-Rpass-analysis=loop-vectorize \
	-framework Accelerate \
	blas2_v1_float.c -o v1_float.o 
	./v1_float.o


# Use accelerate to test
v1a:
	gcc -O3 -march=native -funroll-loops -framework Accelerate -Rpass-missed=loop-vectorize -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize blas2_v1.c -o v1.o 
	./v1.o

v1_:
	gcc-15 -O3 -march=native -ffast-math-fopenmp blas2_v1.c -o v1.o 
	./v1.o

v2:
	gcc -O3 -march=native -funroll-loops -lm \
	-Rpass-missed=loop-vectorize -Rpass=loop-vectorize \
	-Rpass-analysis=loop-vectorize \
	-framework Accelerate \
	blas2_v2.c -o v2.o 
	./v2.o

v2-float:
	gcc -O3 -march=native -funroll-loops -lm \
	-Rpass-missed=loop-vectorize -Rpass=loop-vectorize \
	-Rpass-analysis=loop-vectorize \
	-framework Accelerate \
	blas2_v2_float.c -o v2_float.o 
	./v2_float.o

v2-gcc:
	gcc-15 -O3 -march=native -funroll-loops \
	-fopenmp \
	blas2_v2.c -o v2.o 
	./v2.o