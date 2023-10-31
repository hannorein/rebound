# Miscellaneous tools

## Modulo two pi 
The following function takes the modulo of angle (in radians).
The result is in the range $[0, 2\pi]$.

=== "C"
    ```c
    double f = reb_mod2pi(7.); // returns 0.7168 = 7 - 2*pi
    ```

=== "Python"
    ```python
    f = rebound.mod2pi(7.) // returns 0.7168 = 7 - 2*pi
    ```
## Hash function
REBOUND comes with its own hash function. 
It converts a string to an integer. 
This is used in various parts of the code, mostly to add a more convenient way to refer to particles.
You can call REBOUND's hash function  manually:
=== "C"
    ```c
    uint32_t hash = reb_hash("test string");
    ```

=== "Python"
    ```python
    hash = rebound.hash("test string")
    ```
