Run with `cargo bench -- <bench_name>`.  
Ex: `cargo bench -- qs_32`



Results on Intel(R) Core(TM) i5-10400F CPU @ 2.90GHz, 2904 Mhz, 6 Core(s), 12 Logical Processor(s) running Windows 10, on WSL2:

Single:
```
qs/qs_32                time:   [3.1553 ms 3.1818 ms 3.2086 ms]
qs/qs_64                time:   [7.2810 ms 7.2990 ms 7.3203 ms]
qs/qs_128               time:   [4.0642 s 4.1162 s 4.1746 s]
```

Threaded:
```
qs/qs_threaded_32       time:   [8.2708 ms 8.3255 ms 8.3950 ms]
qs/qs_threaded_64       time:   [10.189 ms 10.551 ms 10.896 ms]
qs/qs_threaded_128      time:   [1.2653 s 1.2798 s 1.2949 s]
qs/qs_threaded_160      time:   [135.69 s 135.89 s 136.07 s]
```

Some hyperparameters (such as sieve_size) or custom bounds / factor bases can be changed to improve the numbers.