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
    if(status == 'PENDING') { message = 'Building...' }
    if(status == 'SUCCESS') { message = 'Build succeeded!' }
    if(status == 'FAILURE') { message = 'Build failed!' }
    if(status == 'ERROR')   { message = 'Build aborted!' }
    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "JenkinsCI/${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: message, state: status]]]])
}

def runCronJob() {
    sh "bash ${HOME}/workspace/build-scripts-cron/cronjob.sh $STAGE_NAME"
}

pipeline {
  agent {
    node { label 'jenkins-slave-1' }
  }
  
  triggers {
    cron('0 3 * * *')
  }
  
  environment {
    SKIP_UPLOAD = '1'
    JOB_TYPE = getJobType()
  }
  
  stages {
    // Stages triggered by GitHub pushes
    stage('notify-pending') {
      when {
        expression { JOB_TYPE == "push" }
      }
      
      steps {
        notifyGitHub('PENDING')
      }
    }
    
    stage('build') {
      when {
        expression { JOB_TYPE == "push" }
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
    
    stage('notify') {
      when {
        expression { JOB_TYPE == "push" }
      }
      
      steps {
        echo 'Setting GitHub status...'
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
          emailext(recipientProviders: [[$class: 'DevelopersRecipientProvider']],  
                  subject: '[JenkinsCI/$PROJECT_NAME] Build # $BUILD_NUMBER - $BUILD_STATUS!', 
                  body: '''${SCRIPT, template="groovy-text.template"}''')
        }
      }
    }
    
    // Stages triggered by cron
    stage('build-scripts-checkout') {
      when {
        expression { JOB_TYPE == "cron" }
      }
      
      steps {
        sh 'cd ${HOME}/workspace/build-scripts-cron/ && git checkout jenkins && git pull --rebase'
      }
    }
    
    stage('centos6') {
      when {
        expression { JOB_TYPE == "cron" }
        expression { SLAVE_OS == "linux" }
      }
      
      steps {
        runCronJob()
      }
    }
    
    stage('centos7') {
      when {
        expression { JOB_TYPE == "cron" }
        expression { SLAVE_OS == "linux" }
      }
      
      steps {
        runCronJob()
      }
    }
    
    stage('mac') {
      when {
        expression { JOB_TYPE == "cron" }
        expression { SLAVE_OS == "mac" }
      }
      
      steps {
        runCronJob()
      }
    }
    
    stage('build-scripts-reset') {
      when {
        expression { JOB_TYPE == "cron" }
      }
      
      steps {
        sh 'cd ${HOME}/workspace/build-scripts-cron/ && git checkout master'
      }
    }
  }
}
