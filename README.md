
### How to compile

```
. env.sh  ## load module
cd ./fastF
export CC=/path/to/gcc
cmake -B ./build -S .
```

### How to ignore one folder or file already in the git repository

1. add the `folder_name` or `file_name` into `.gitignore`
2. `git rm -r --cached <folder/file_name>`




