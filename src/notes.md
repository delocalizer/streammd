## Description

As this is partly a learning exercise I'm going to drop some notes and
references that I found useful about style, idiom etc.

## Unwinding
https://doc.rust-lang.org/nomicon/unwinding.html
> Rust has a tiered error-handling scheme:
> * If something might reasonably be absent, Option is used.
> * If something goes wrong and can reasonably be handled, Result is used.
> * If something goes wrong and cannot reasonably be handled, the thread panics.
> * If something catastrophic happens, the program aborts.
>
> ...
>
> Don't build your programs to unwind under normal circumstances. Ideally, you should only panic for programming errors or extreme problems.


## Mapping modules to the filesystem

http://web.mit.edu/rust-lang_v1.25/arch/amd64_ubuntu1404/share/doc/rust/html/book/second-edition/ch07-01-mod-and-the-filesystem.html

In particular [moving modules to other files](http://web.mit.edu/rust-lang_v1.25/arch/amd64_ubuntu1404/share/doc/rust/html/book/second-edition/ch07-01-mod-and-the-filesystem.html#moving-modules-to-other-files)

## Idiomatic error handling

* Use [?](https://doc.rust-lang.org/reference/expressions/operator-expr.html#the-question-mark-operator)
  wherever possible on `Result` or `Option` types. This expands to `expr`
  for the success branch (`Ok(expr)` | `Some(expr)`), and an early return of
  `Err(err)` | `None` from the 'failure' branch. This is much more elegant than
  simply `unwrap`ing `Result` or `Option` to get a value, which will `panic` on
  the failure branch. It is also generally more elegant than chaining together 
  multiple `match` expressions to propagate `Result`s or `Option`s through a
  function to get the final return type.
* Using the [anyhow](https://docs.rs/anyhow/latest/anyhow/) crate simplifies
  error handling. You can have a generic `Result<T>` return type instead of a
  specific `Result<T, E>` or clunky `Result<T, Box<dyn std::error:Error>>`
  for functions that can return multiple error types. It also provides the
  `.context()` and `.with_context` methods to provide information at the point
  of the error.
