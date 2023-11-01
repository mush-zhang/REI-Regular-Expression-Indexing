#!/bin/bash

num_ngrams_list=( 4 8 16 32 64 96 128 192 256 320 384 448 512 )
grams_list=( bigram )

rm index_building.csv
touch index_building.csv

for num_ngram in ${num_ngrams_list[*]}; do
    for grams in ${grams_list[*]}; do
        rm -r new_${grams}_${num_ngram}
        rm ../microbench${num_ngram}.o
        g++ -O3 -std=c++17 -Ofast -march=native -mfma -mavx -fomit-frame-pointer -ffp-contract=fast -flto -DARMA_NO_DEBUG -pthread ${grams}_record_all_${num_ngram}.cpp  -L/usr/local/lib/ -lre2 -lstdc++fs -o ../${grams}_record_all_${num_ngram}.o
    done
done 
cd ..

for num_ngram in ${num_ngrams_list[*]}; do
    for grams in ${grams_list[*]}; do
        ./${grams}_record_all_${num_ngram}.o 0 
    done
done

cd microbench_script
