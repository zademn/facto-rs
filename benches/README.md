Run with `cargo bench -- <bench_name>`.  
Ex: `cargo bench -- qs_32`



Results on Intel(R) Core(TM) i5-10400F CPU @ 2.90GHz, 2904 Mhz, 6 Core(s), 12 Logical Processor(s) running Windows 10, on WSL2:

Single:
```
qs/qs_32                time:   [4.2872 ms 4.2967 ms 4.3134 ms]
qs/qs_64                time:   [9.9433 ms 9.9932 ms 10.036 ms]
qs/qs_128               time:   [4.4153 s 4.4337 s 4.4565 s]
```

Threaded:
```
qs/qs_threaded_32       time:   [11.147 ms 11.428 ms 11.588 ms]
qs/qs_threaded_64       time:   [18.051 ms 18.841 ms 20.146 ms]
qs/qs_threaded_128      time:   [1.5326 s 1.5560 s 1.5862 s]
qs/qs_threaded_160      time:   [154.85 s 156.94 s 159.55 s]
```

Some hyperparameters (such as sieve_size) or custom bounds / factor bases can be changed to improve the numbers.