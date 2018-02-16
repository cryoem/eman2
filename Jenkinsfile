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
    if(JOB_TYPE == "push") {
        if(status == 'PENDING') { message = 'Building...' }
        if(status == 'SUCCESS') { message = 'Build succeeded!' }
        if(status == 'FAILURE') { message = 'Build failed!' }
        if(status == 'ERROR')   { message = 'Build aborted!' }
        step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "JenkinsCI/${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: message, state: status]]]])
    }
}

def notifyEmail() {
    if(JOB_TYPE == "push") {
        emailext(recipientProviders: [[$class: 'DevelopersRecipientProvider']],  
                 subject: '[JenkinsCI/$PROJECT_NAME/push] ' + "($GIT_BRANCH_SHORT - ${GIT_COMMIT_SHORT})" + ' #$BUILD_NUMBER - $BUILD_STATUS!',
                 body: '''${SCRIPT, template="groovy-text.template"}''',
                 attachLog: true
                 )
    }
}

def isRelease() {
    return GIT_BRANCH ==~ /.*\/release.*/
}

def isContinuousBuild() {
    return GIT_BRANCH ==~ /.*\/master/
}

def stage_name_to_os(stage_name) {
    def result = ['centos6': 'linux',
                  'centos7': 'linux',
                  'mac':     'mac',
                  'win':     'win'
                  ]
    
    return result[stage_name]
}

def isBuildStage() {
    def buildStage = ['centos6': CI_BUILD_LINUX,
                      'centos7': CI_BUILD_LINUX,
                      'mac':     CI_BUILD_MAC,
                      'win':     CI_BUILD_WIN
                      ]
    
    return (stage_name_to_os(STAGE_NAME) == SLAVE_OS && (CI_BUILD == "1" || buildStage[STAGE_NAME] == "1"))
}

def isRunCurrentStage() {
    return isContinuousBuild() && isBuildStage()
}

def runCronJob() {
    if(isContinuousBuild()) {
        if(SLAVE_OS != 'win')
            sh "rsync -avzh --stats ${INSTALLERS_DIR}/eman2.${SLAVE_OS}.sh ${DEPLOY_DEST}/eman2.${STAGE_NAME}.unstable.sh"
        else
            bat 'ci_support\\rsync_wrapper.bat'
    }
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
    node { label 'jenkins-slave-1' }
  }
  
  options {
    disableConcurrentBuilds()
    timestamps()
  }
  
  environment {
    JOB_TYPE = getJobType()
    GIT_BRANCH_SHORT = sh(returnStdout: true, script: 'echo ${GIT_BRANCH##origin/}').trim()
    GIT_COMMIT_SHORT = sh(returnStdout: true, script: 'echo ${GIT_COMMIT:0:7}').trim()
    HOME_DIR = getHomeDir()
    INSTALLERS_DIR = '${HOME_DIR}/workspace/${STAGE_NAME}-installers'
    DEPLOY_DEST    = 'zope@ncmi.grid.bcm.edu:/home/zope/zope-server/extdata/reposit/ncmi/software/counter_222/software_136/'
    NUMPY_VERSION='1.9'

    CI_BUILD       = sh(script: "! git log -1 | grep '.*\\[ci build\\].*'",       returnStatus: true)
    CI_BUILD_WIN   = sh(script: "! git log -1 | grep '.*\\[ci build win\\].*'",   returnStatus: true)
    CI_BUILD_LINUX = sh(script: "! git log -1 | grep '.*\\[ci build linux\\].*'", returnStatus: true)
    CI_BUILD_MAC   = sh(script: "! git log -1 | grep '.*\\[ci build mac\\].*'",   returnStatus: true)
  }
  
  stages {
    // Stages triggered by GitHub pushes
    stage('notify-pending') {
      steps {
        notifyGitHub('PENDING')
      }
    }
    
    stage('build') {
      when {
        not { expression { isContinuousBuild() } }
        expression { false }
      }
      
      parallel {
        stage('recipe') {
          steps {
            sh 'bash ci_support/build_recipe.sh'
          }
        }
        
        stage('no_recipe') {
          steps {
            sh 'source $(conda info --root)/bin/activate eman-env && bash ci_support/build_no_recipe.sh'
          }
        }
      }
    }
    
    stage('centos6') {
      when {
        expression { isRunCurrentStage() }
      }
      
      steps {
        runCronJob()
      }
    }
    
    stage('centos7') {
      when {
        expression { isRunCurrentStage() }
      }
      
      steps {
        runCronJob()
      }
    }
    
    stage('mac') {
      when {
        expression { isRunCurrentStage() }
      }
      environment {
        EMAN_TEST_SKIP=1
      }
      
      steps {
        runCronJob()
      }
    }
    
    stage('win') {
      when {
        expression { isRunCurrentStage() }
      }
      
      steps {
        runCronJob()
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
