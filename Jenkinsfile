pipeline {
  agent {
    node {
      label 'jenkins-slave-1'
    }
  }
  
  environment {
    SKIP_UPLOAD = '1'
  }
  
  stages {
    stage('recipe') {
      steps {
        sh 'bash ci_support/build_recipe.sh'
      }
    }
    
    stage('no_recipe') {
      steps {
        sh 'bash ci_support/build_no_recipe.sh'
      }
    }
  }
}
