# 0) Packages
install.packages(c("usethis", "gitcreds"))
library(usethis)
library(gitcreds)

# 1) Make sure Git exists
system("git --version")

# 2) Set your Git identity (one-time global)
use_git_config(user.name = "isaactklein", user.email = "isaactklein@gmail.com")

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                # 3) Initialize Git in this project (creates .git, first commit)
use_git()   # RStudio will ask to restart; say yes.

# 4) Create a GitHub token and store it (one-time)
create_github_token()      # opens browser → generate token with "repo" scope
gitcreds::gitcreds_set()             # paste token here when prompted

# 5) Push this project to GitHub (creates repo + adds 'origin' + pushes main)
use_github(private = FALSE) # or FALSE if you want public
system("git --version")   # should print a version
Sys.which("git")          # should show a path like "C:\\Program Files\\Git\\cmd\\git.exe"
usethis::git_sitrep()     # optional, good diagnostics

# Try common install locations and prepend to PATH if found
cands <- c("C:/Program Files/Git/cmd",
           "C:/Program Files/Git/bin",
           "C:/Program Files (x86)/Git/cmd")
hit <- cands[dir.exists(cands)]
if (length(hit)) Sys.setenv(PATH = paste(hit[1], Sys.getenv("PATH"), sep=";"))

system('git remote -v')         # should show https://github.com/isaactklein/capstone_project.git
system('git status')            # should say "nothing to commit" before you push
system('git branch -vv')        # see which branch you’re on (you’re on master)

system('git push -u origin master')   # push your "Initial commit"


# 1) Set identity for THIS repository
usethis::use_git_config(
  user.name = "isaactklein",
  user.email = "isaactklein@gmail.com",
  scope = "project"
)

# sanity check
system('git config user.name'); system('git config user.email')

# 2) Commit whatever’s pending
system('git add -A')
system('git commit -m "Sync to friend repo"')

# 3) Push your current branch to your friend's repo on branch 'isaac'
usethis::use_git_remote("friend", "https://github.com/JohnPaulCAcosta/Drug-Repurposing-Capstone.git", overwrite = TRUE)
br <- system('git branch --show-current', intern = TRUE)
system(paste("git push -u friend", br))

# 4) Open PR page
browseURL("https://github.com/JohnPaulCAcosta/Drug-Repurposing-Capstone/pull/new/isaac")


  
  usethis::create_github_token()
gitcreds::gitcreds_set()


