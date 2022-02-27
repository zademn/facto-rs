Run with `cargo bench -- <bench_name>`.  
Ex: `cargo bench -- qs_32`



Results on Intel(R) Core(TM) i5-10400F CPU @ 2.90GHz, 2904 Mhz, 6 Core(s), 12 Logical Processor(s) running Manjaro KDE:

Single:
```
single/qs_32            time:   [2.9331 ms 2.9614 ms 3.0210 ms]
single/qs_64            time:   [4.9804 ms 5.0204 ms 5.0529 ms]
single/qs_128           time:   [2.2731 s 2.2870 s 2.3029 s]
```

Threaded:
```
qs/qs_threaded_32       time:   [8.4680 ms 8.8760 ms 9.3579 ms]
qs/qs_threaded_64       time:   [9.8993 ms 10.133 ms 10.630 ms]
qs/qs_threaded_128      time:   [837.28 ms 843.95 ms 850.63 ms]
qs/qs_threaded_160      time:   [127.83 s 128.41 s 129.34 s]
```

Some hyperparameters (such as sieve_size) or custom bounds / factor bases can be changed to improve the numbers.