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
# 
# model.schizo2$importance
# 0            1 MeanDecreaseAccuracy
# xlogp                         0.0018641831  0.017919565          0.003553538
# tpsa                          0.0043039582  0.011116356          0.005048575
# num.atoms                     0.0043610876  0.054934154          0.010313051
# mass                          0.0051332228  0.100722032          0.016317274
# DRD                           0.0076452204  0.150752700          0.024299748
# HTR                           0.0027154958 -0.008158837          0.001404785
# HRH                           0.0001204028  0.019893134          0.002333436
# dopamine receptor antagonist  0.0241921568  0.328116707          0.060014962
# serotonin receptor antagonist 0.0001839752  0.011689649          0.001567902
# DRD.count                     0.0102778172  0.142230392          0.025914419
# HTR.count                     0.0032516486  0.006212620          0.003563729
# fp_151                        0.0008337016  0.005290095          0.001274531
# MeanDecreaseGini
# xlogp                                4.9936445
# tpsa                                 3.0194489
# num.atoms                            3.2713608
# mass                                 5.0114125
# DRD                                  7.2148980
# HTR                                  1.7655458
# HRH                                  1.1128437
# dopamine receptor antagonist        16.8630613
# serotonin receptor antagonist        1.2540337
# DRD.count                            7.6862336
# HTR.count                            2.6632006
# fp_151                               0.7477393
# > schizo.test2 = model.schizo2$test
# > confusionMatrix(schizo.test2$predicted, as.factor(testing.drugs$schizophrenia), positive = "1")
# Confusion Matrix and Statistics
# 
# Reference
# Prediction  0  1
# 0 99  3
# 1  3 11
# 
# Accuracy : 0.9483          
# 95% CI : (0.8908, 0.9808)
# No Information Rate : 0.8793          
# P-Value [Acc > NIR] : 0.01027         
# 
# Kappa : 0.7563          
# 
# Mcnemar's Test P-Value : 1.00000         
#                                           
#             Sensitivity : 0.78571         
#             Specificity : 0.97059         
#          Pos Pred Value : 0.78571         
#          Neg Pred Value : 0.97059         
#              Prevalence : 0.12069         
#          Detection Rate : 0.09483         
#    Detection Prevalence : 0.12069         
#       Balanced Accuracy : 0.87815         
#                                           
#        'Positive' Class : 1               
#                                           
# > confusionMatrix(model.schizo2$predicted, model.schizo2$y, positive = "1")
# Confusion Matrix and Statistics
# 
#           Reference
# Prediction   0   1
#          0 231   7
#          1   5  25
#                                           
#                Accuracy : 0.9552          
#                  95% CI : (0.9231, 0.9767)
#     No Information Rate : 0.8806          
#     P-Value [Acc > NIR] : 2.042e-05       
#                                           
#                   Kappa : 0.7812          
#                                           
#  Mcnemar's Test P-Value : 0.7728          
# 
# Sensitivity : 0.78125         
# Specificity : 0.97881         
# Pos Pred Value : 0.83333         
# Neg Pred Value : 0.97059         
# Prevalence : 0.11940         
# Detection Rate : 0.09328         
# Detection Prevalence : 0.11194         
# Balanced Accuracy : 0.88003         
# 
# 'Positive' Class : 1               
# 
# > # 2) Commit whatever’s pending
