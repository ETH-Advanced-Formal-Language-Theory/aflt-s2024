# Advanced Formal Language Theory (263-5352-00L; Spring 2024)

[Website](https://rycolab.io/classes/aflt-s24/)

This course serves as an introduction to various advanced topics in formal language theory. The primary focus of the course is on weighted formalisms, which can easily be applied in machine learning. Topics include finite-state machines as well as the algorithms that are commonly used for their manipulation. We will also cover weighted context-free grammars, weighted tree automata, and weighted mildly context-sensitive formalisms.

## About the Repository
This repository will contain pieces of the `rayuela` library and some teaching material which you'll need to complete the course assignments.

## Installation
To install the course library run 
```
bash autograding_tests/setup_test_env.sh
```
in the base directory.

## Code Updates
We are constantly working to improve `rayuela`! Unfortunately, this means that sometimes we will have to update the skeleton structure in the template repository.
There are two possibilities to merge your private repository with the updated template:
- re-accepting the assignments;
- merging the remote template code.

In both cases, you won't be able to automatically merge the template code with yours, as your implementations will always create merge conflict.
Complete the merging with care, as you might lose part of your code or important updates to the template!

#### Re-accepting the Assignments
This can be done by deleting your private repository (on Github) and joining Classroom again [here](https://classroom.github.com/a/_YxlmtX1).
After you re-accepted the assignment, you will have to locally pull and merge the updated template with your repo with
```
git pull --allow-unrelated-histories
```

#### Merging the Template
To merge the latest changes in the template repository to your private repository, use
```
git remote add template https://github.com/ETH-Advanced-Formal-Language-Theory/aflt-s2024
git fetch --all
git merge template/main --allow-unrelated-histories
```


## Testing
Every commit that you push to this repository will be automatically evaluated by a dedicated Github Workflow.
For debugging your implementations, you can also run the tests locally. To run all the tests, you can use `pytest .` in the base directory. If you want to test only a specific assignment, you can use `pytest autograding_tests/test_hw[num].py`.
