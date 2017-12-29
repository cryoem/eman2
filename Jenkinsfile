def getJobType() {
    def causes = "${currentBuild.rawBuild.getCauses()}"
    def job_type = "UNKNOWN"
    
    if(causes ==~ /.*TimerTrigger.*/) {
        job_type = "cron"
    }
    if(causes ==~ /.*GitHubPushCause.*/) {
        job_type = "push"
    }
    if(causes ==~ /.*UserIdCause.*/) {
        job_type = "manual"
    }
    if(causes ==~ /.*ReplayCause.*/) {
        job_type = "manual"
    }
    
    return job_type
}

def notifyGitHub(status) {
    if(status == 'PENDING') {
        message = 'Building...'
    }
    if(status == 'SUCCESS') {
        message = 'Build succeeded!'
    }
    if(status == 'FAILURE') {
        message = 'Build failed!'
    }
    if(status == 'ERROR') {
        message = 'Build aborted!'
    }
    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: message, state: status]]]])
}

pipeline {
  agent {
    node {
      label 'jenkins-slave-1'
    }
  }
  
  environment {
    SKIP_UPLOAD = '1'
    JOB_TYPE = getJobType()
  }
  
  stages {
    stage('pending') {
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
      }
    }
  }
}
