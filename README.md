# REI-Regular-Expression-Indexing
REI: Regular Expression Indexing

REI is a framework for indexing regular expressions. We use an n-gram-based indexing strategy and an efficient storage mechanism that results in a speedup of up to 14x compared to state-of-the-art regex processing engines that do not use indexing, using only 2.1% of extra space. REI is also modular and can work with existing regular expression packages, making it easy to deploy in a variety of settings.

**This project is under active development**

## Prerequisites

We use Google-RE2 as for full regex matching after index filtering.

- Google-RE2

    Please follow the instruction in [RE2 Github Repository](https://github.com/google/re2).

REI is implemented in C++, and we provide cmake file for building the project with the external dependencies excluding `g++` and `Boost`. For Ubuntu as an example, you can do

```bash
sudo apt update && sudo apt upgrade
```

- g++ (version 8.4.0 or higher)
  
    ```bash
    sudo apt install build-essential
    ```

- Boost Library (version 1.65.1.0 or higher)
  
    ```bash
    sudo apt-get install libboost-all-dev
    ```

Make sure to check if the version satifies the requirement by

```bash
g++ --version
cmake --version
dpkg -s libboost-all-dev | grep 'Version'
```

US-Accident dataset is accessible [**here**](https://github.com/mush-zhang/Blare/tree/main/data)

## Experiments

1. How does the chosen value of 𝑞 for 𝑛-grams of the index affect the overhead of constructing the index and the improvement of performance? [**In ngram_w_inverted**](https://github.com/mush-zhang/REI-Regular-Expression-Indexing/blob/main/ngram_w_inverted/ngram_type.ipynb)
2. How does the chosen value of 𝑘 (the number of 𝑛-grams used in the index) affect the overhead of constructing the index and the improvement of performance? [**In ngram_w_inverted**](https://github.com/mush-zhang/REI-Regular-Expression-Indexing/blob/main/ngram_w_inverted/bigram_num.ipynb)
3. How does the existence of prior knowledge of regex workload impact the performance improvement of the constructed index? [**In ngram_english**](https://github.com/mush-zhang/REI-Regular-Expression-Indexing/tree/main/ngram_english)
4. How does the granularity of the index affect the overhead of constructing the index and the improvement of performance? [**In ngram_multilevel_inverted**](https://github.com/mush-zhang/REI-Regular-Expression-Indexing/blob/main/ngram_multilevel_inverted/multilevel_bitvect_only.ipynb)
5. How does REI compare with other commonly used indexing schema on the overhead of constructing the index and the improvement of performance? [**In ngram_multilevel_inverted**](https://github.com/mush-zhang/REI-Regular-Expression-Indexing/blob/main/ngram_multilevel_inverted/multilevel_128only.ipynb)