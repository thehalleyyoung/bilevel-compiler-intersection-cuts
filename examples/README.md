# Examples

`simple_bilevel.toml` is a bilevel LP with one leader variable and one
follower variable:

```
min_x  -x - 7y
s.t.   -2x + y <= 4,  0 <= x <= 10
       y in argmin_y { -y : -x + y <= 1,  x + y <= 5,  0 <= y <= 10 }
```

Its bilevel optimum is x = 2, y = 3 with leader objective -23.

```bash
bicut compile --input examples/simple_bilevel.toml --output simple.mps
bicut compile --input examples/simple_bilevel.toml --output simple.lp
```

Both commands write the high-point relaxation (all constraints, follower
optimality dropped). For this instance its optimum coincides with the
bilevel optimum; HiGHS returns -23 at (2, 3) on either file.
