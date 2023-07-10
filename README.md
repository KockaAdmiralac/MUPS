# MUPS
Homework for our [Multiprocessor Systems class](http://mups.etf.bg.ac.rs/) at [School of Electrical Engineering, University of Belgrade](https://www.etf.bg.ac.rs/). Worked on by [@topofkeks](https://github.com/topofkeks) and I, we had three homework tasks that involved programming using [OpenMP](https://en.wikipedia.org/wiki/OpenMP), [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) and [CUDA](https://en.wikipedia.org/wiki/CUDA), and one that involved a [cache coherency simulator](http://mups.etf.bg.ac.rs/simulatori/). The programming ones involved programs that needed to be sped up using the mentioned technologies. The simulator one was individual, and the report was supplied in a text file.

Our reports can be found in [Releases](https://github.com/KockaAdmiralac/MUPS/releases).

## Homework
All homework tasks are in Serbian.

1. [First homework task (OpenMP)](https://web.archive.org/web/20230710231639im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ1_2022-2023.pdf)
    - [Initial code for the first task](https://web.archive.org/web/20230710231703im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ1_OpenMP.zip)
2. [Second homework task (MPI)](https://web.archive.org/web/20230710231655im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ2_2022-2023.pdf)
    - [Initial code for the second task](https://web.archive.org/web/20230710231737im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ2_MPI.zip)
3. [Third homework task (cache coherency simulation)](https://web.archive.org/web/20230710231651im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ3_2022-2023.pdf)
4. [Fourth homework task (CUDA)](https://web.archive.org/web/20230710231726im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ4_2022-2023.pdf)
    - [Initial code for the fourth task](https://web.archive.org/web/20230710231717im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ4_CUDA.zip)

## Structure
- `run.py` is our test runner that plots our speedup graphs and logs the output
- `commit.sh` is a script that auto-commits to our SVN repositories (the homework had to go through SVN for some reason). For that, it needed the following environment variables:
    - `SSH_HOST`: Host to connect to through SSH
    - `SSH_PORT`: SSH port
    - `SSH_USERNAME`: Username to use when connecting through SSH
    - `SSH_KEY`: SSH password
    - `REPO_LOCATION`: Our Git repository location
- `Makefile` builds all tasks for us
- `src` contains homework files named as `dzXzY.c[u]`, where X is the homework task number and Y is the subtask number
