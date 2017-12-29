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
      steps {
        step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'Building...', state: 'PENDING']]]])
      }
    }
    
    stage('build') {
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
  }
  
  post {
    success {
      step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'Build succeeded!', state: 'SUCCESS']]]])
    }
    
    failure {
      step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'Build failed!', state: 'FAILURE']]]])
    }
    
    aborted {
      step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "${JOB_NAME}"], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'Build aborted!', state: 'ERROR']]]])
    }
  }
}
