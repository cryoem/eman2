def getJobType() {
    def causes = "${currentBuild.rawBuild.getCauses()}"
    def job_type = "UNKNOWN"
    
    if(causes ==~ /.*TimerTrigger.*/)    { job_type = "cron" }
    if(causes ==~ /.*GitHubPushCause.*/) { job_type = "push" }
    if(causes ==~ /.*UserIdCause.*/)     { job_type = "manual" }
    if(causes ==~ /.*ReplayCause.*/)     { job_type = "manual" }
    
    return job_type
}

def notifyGitHub(status) {
    if(JOB_TYPE == "push" || NOTIFY_GITHUB == "true") {
        if(status == 'PENDING') { message = 'Stage: ' + (env.PARENT_STAGE_NAME ?: STAGE_NAME) }
        if(status == 'SUCCESS') { message = 'Build succeeded!' }
        if(status == 'FAILURE') { message = 'Build failed!' }
        if(status == 'ERROR')   { message = 'Build aborted!' }
        step([$class: 'GitHubCommitStatusSetter', 
              contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "JenkinsCI/${JOB_NAME}"], 
              statusResultSource: [$class: 'ConditionalStatusResultSource', 
                                   results: [[$class: 'AnyBuildResult', message: message, state: status]]]])
    }
}

def notifyEmail() {
    if(JOB_TYPE == "push" || NOTIFY_EMAIL == "true") {
        emailext(to: "$GIT_AUTHOR_EMAIL",  
                 subject: '[JenkinsCI/$PROJECT_NAME] ' + "($GIT_BRANCH_SHORT - ${GIT_COMMIT_SHORT})" + ' #$BUILD_NUMBER - $BUILD_STATUS!',
                 body: '''${SCRIPT, template="groovy-text.template"}''',
                 attachLog: true
                 )
    }
}

def selectNotifications() {
    if(env.JOB_TYPE == 'manual') {
        def result = input(message: 'Select notifications:',
                           parameters :
                                   [booleanParam(defaultValue: false, description: 'Notify GitHub?', name: 'notify_github'),
                                    booleanParam(defaultValue: false, description: 'Email author?',  name: 'notify_email')]
                           )
                 
        env.NOTIFY_GITHUB = result.notify_github
        env.NOTIFY_EMAIL  = result.notify_email
    }
    else {
        env.NOTIFY_GITHUB = true
        env.NOTIFY_EMAIL  = true
    }
}

def isMasterBranch() {
    return GIT_BRANCH_SHORT == "master"
}

def isReleaseBranch() {
    return GIT_BRANCH_SHORT ==~ /release.*/
}

def isContinuousBuild() {
    return (CI_BUILD == "1" && isMasterBranch()) || isReleaseBranch()
}

def isExperimentalBuild() {
    return CI_BUILD == "1" && !(isMasterBranch() || isReleaseBranch())
}

def isBinaryBuild() {
    return isContinuousBuild() || isExperimentalBuild()
}

def testPackage(suffix, dir) {
    if(SLAVE_OS != 'win')
        sh "bash tests/test_binary_installation.sh ${INSTALLERS_DIR}/eman2" + suffix + ".${SLAVE_OS}.sh ${INSTALLERS_DIR}/" + dir
    else
        sh 'ci_support/test_wrapper.sh ' + suffix + ' ' + dir
}

def deployPackage() {
    if(isContinuousBuild()) {
        upload_dir = 'continuous_build'
        upload_ext = 'unstable'
    }
    if(isExperimentalBuild()) {
        upload_dir = 'experimental'
        upload_ext = 'experimental'
    }
    
    if(SLAVE_OS != 'win') {
        sh "rsync -avzh --stats ${INSTALLERS_DIR}/eman2.${SLAVE_OS}.sh ${DEPLOY_DEST}/" + upload_dir + "/eman2." + JOB_NAME.toLowerCase() + "." + upload_ext + ".sh"
        sh "rsync -avzh --stats ${INSTALLERS_DIR}/eman2_huge.${SLAVE_OS}.sh ${DEPLOY_DEST}/" + upload_dir + "/eman2_huge." + JOB_NAME.toLowerCase() + "." + upload_ext + ".sh"
    }
    else
        bat 'ci_support\\rsync_wrapper.bat ' + upload_dir + ' ' + upload_ext
}

def getHomeDir() {
    def result = ''
    if(SLAVE_OS == "win") {
        result = "${USERPROFILE}"
    }
    else {
        result = "${HOME}"
    }
    
    return result
}

pipeline {
  agent {
    node { label "${JOB_NAME}-slave" }
  }
  
  options {
    disableConcurrentBuilds()
    timestamps()
  }
  
  environment {
    JOB_TYPE = getJobType()
    GIT_BRANCH_SHORT = sh(returnStdout: true, script: 'echo ${GIT_BRANCH##origin/}').trim()
    GIT_COMMIT_SHORT = sh(returnStdout: true, script: 'echo ${GIT_COMMIT:0:7}').trim()
    GIT_AUTHOR_EMAIL = sh(returnStdout: true, script: 'git log -1 --format="%ae"').trim()
    HOME_DIR = getHomeDir()
    HOME = "${HOME_DIR}"     // on Windows HOME is set to something like C:\Program Files\home\eman
    INSTALLERS_DIR = '${HOME_DIR}/workspace/${JOB_NAME}-installers'

    CI_BUILD       = sh(script: "! git log -1 | grep '.*\\[ci build\\].*'",       returnStatus: true)
  }
  
  stages {
    stage('init') {
      options {
        timeout(time: 10, unit: 'MINUTES') 
      }
      
      steps {
        selectNotifications()
        notifyGitHub('PENDING')
        sh 'env | sort'
      }
    }
    
    stage('build-local') {
      when {
        not { expression { isBinaryBuild() } }
        expression { JOB_NAME != 'Win' }
      }
      
      steps {
        notifyGitHub('PENDING')
        sh 'source $(conda info --root)/bin/activate eman-deps-14.0 && bash ci_support/build_no_recipe.sh'
      }
    }
    
    stage('build-recipe') {
      steps {
        notifyGitHub('PENDING')
        sh 'bash ci_support/build_recipe.sh'
      }
    }
    
    stage('package') {
      when {
        expression { isBinaryBuild() }
      }
      environment {
        PARENT_STAGE_NAME = "${STAGE_NAME}"
      }
      
      parallel {
        stage('notify') {
          steps {
            notifyGitHub('PENDING')
          }
        }
        stage('mini') {
          steps {
            sh "bash ci_support/package.sh ${INSTALLERS_DIR} " + '${WORKSPACE}/ci_support/constructor-mini/'
          }
        }
        stage('huge') {
          steps {
            sh "bash ci_support/package.sh ${INSTALLERS_DIR} " + '${WORKSPACE}/ci_support/constructor-huge/'
          }
        }
      }
    }
    
    stage('test-package') {
      when {
        expression {isBinaryBuild() }
      }
      environment {
        PARENT_STAGE_NAME = "${STAGE_NAME}"
      }
      
      parallel {
        stage('notify') {
          steps {
            notifyGitHub('PENDING')
          }
        }
        stage('mini') {
          steps {
            testPackage('', 'mini')
          }
        }
        stage('huge') {
          steps {
            testPackage('_huge','huge')
          }
        }
      }
    }
    
    stage('deploy') {
      when {
        expression {isBinaryBuild() }
      }
      
      steps {
        notifyGitHub('PENDING')
        deployPackage()
      }
    }
  }
  
  post {
    success {
      notifyGitHub('SUCCESS')
    }
    
    failure {
      notifyGitHub('FAILURE')
    }
    
    aborted {
      notifyGitHub('ERROR')
    }
    
    always {
      notifyEmail()
    }
  }
}
