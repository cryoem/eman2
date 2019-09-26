def binary_size_suffix = ['mini':'', 'huge':'_huge']

def getOSName() {
    if(!isUnix()) return 'win'
    else {
        uname  = sh(returnStdout: true, script: 'uname -s').trim().toLowerCase()
        os_map = ['linux':'linux', 'darwin':'mac']

        return os_map[uname]
    }
}

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
        if(status == 'ABORTED') { message = 'Build aborted!'; status == 'ERROR' }
        step([$class: 'GitHubCommitStatusSetter', 
              contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: "JenkinsCI/${AGENT_OS_NAME.capitalize()}"],
              statusResultSource: [$class: 'ConditionalStatusResultSource', 
                                   results: [[$class: 'AnyBuildResult', message: message, state: status]]]])
    }
}

def notifyEmail() {
    if(JOB_TYPE == "push" || NOTIFY_EMAIL == "true") {
        emailext(to: "$GIT_AUTHOR_EMAIL",
                 from: "JenkinsCI ($AGENT_OS_NAME) <jenkins@jenkins>",
                 subject: '$BUILD_STATUS! ' + "($GIT_BRANCH_SHORT - ${GIT_COMMIT_SHORT})" + ' #$BUILD_NUMBER',
                 body: '''${SCRIPT, template="groovy-text.template"}''',
                 attachLog: true
                 )
    }
    if(JOB_TYPE == "cron") {
        emailext(to: '$DEFAULT_RECIPIENTS',
                 from: "JenkinsCI ($SLAVE_OS) <jenkins@jenkins>",
                 subject: '[cron] - $BUILD_STATUS! ' + "($GIT_BRANCH_SHORT - ${GIT_COMMIT_SHORT})" + ' #$BUILD_NUMBER',
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
    else if(env.JOB_TYPE == 'cron') {
        env.NOTIFY_GITHUB = false
        env.NOTIFY_EMAIL  = false
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
    return (CI_BUILD == "1" && isMasterBranch()) || isReleaseBranch() || JOB_TYPE == "cron"
}

def isExperimentalBuild() {
    return CI_BUILD == "1" && !(isMasterBranch() || isReleaseBranch())
}

def isBinaryBuild() {
    return isContinuousBuild() || isExperimentalBuild()
}

def testPackage(suffix, dir) {
    if(isUnix())
        sh  "bash tests/test_binary_installation.sh   ${WORKSPACE}/eman2"  + suffix + ".${AGENT_OS_NAME}.sh ${INSTALLERS_DIR}/"  + dir
    else
        bat "call tests\\test_binary_installation.bat ${WORKSPACE}\\eman2" + suffix + ".win.exe        ${INSTALLERS_DIR}\\" + dir
}

def deployPackage(size_type='') {
    if(isContinuousBuild())   stability_type = 'unstable'
    if(isExperimentalBuild()) stability_type = 'experimental'

    if(isUnix()) installer_ext = 'sh'
    else         installer_ext = 'exe'

    sshPublisher(publishers: [
                              sshPublisherDesc(configName: 'Installer-Server',
                                               transfers:
                                                          [sshTransfer(sourceFiles: "eman2" + size_type + ".${AGENT_OS_NAME}." + installer_ext,
                                                                       removePrefix: "",
                                                                       remoteDirectory: stability_type,
                                                                       remoteDirectorySDF: false,
                                                                       cleanRemote: false,
                                                                       excludes: '',
                                                                       execCommand: "cd ${DEPLOY_PATH}/" + stability_type + " && mv eman2" + size_type + ".${AGENT_OS_NAME}." + installer_ext + " eman2" + size_type + ".${AGENT_OS_NAME}." + stability_type + "." + installer_ext,
                                                                       execTimeout: 120000,
                                                                       flatten: false,
                                                                       makeEmptyDirs: false,
                                                                       noDefaultExcludes: false,
                                                                       patternSeparator: '[, ]+'
                                                                      )
                                                          ],
                                                          usePromotionTimestamp: false,
                                                          useWorkspaceInPromotion: false,
                                                          verbose: true
                                              )
                             ]
                )
}

def getHomeDir() {
    if(!isUnix()) return "${USERPROFILE}"
    else          return "${HOME}"
}

pipeline {
  agent {
    node { label "eman-build-agent" }
  }
  
  options {
    disableConcurrentBuilds()
    timestamps()
  }
  
  environment {
    AGENT_OS_NAME = getOSName()
    JOB_TYPE = getJobType()
    GIT_BRANCH_SHORT = sh(returnStdout: true, script: 'echo ${GIT_BRANCH##origin/}').trim()
    GIT_COMMIT_SHORT = sh(returnStdout: true, script: 'echo ${GIT_COMMIT:0:7}').trim()
    GIT_AUTHOR_EMAIL = sh(returnStdout: true, script: 'git log -1 --format="%ae"').trim()
    HOME_DIR = getHomeDir()
    HOME = "${HOME_DIR}"     // on Windows HOME is set to something like C:\Program Files\home\eman
    INSTALLERS_DIR = sh(returnStdout: true, script: "python -c 'import os; print(os.path.join(\"${HOME_DIR}\", \"workspace\", \"jenkins-eman-installers\"))'").trim()

    CI_BUILD       = sh(script: "! git log -1 | grep '.*\\[ci build\\].*'",       returnStatus: true)
  }
  
  stages {
    stage('init') {
      options { timeout(time: 10, unit: 'MINUTES') }
      
      steps {
        selectNotifications()
        notifyGitHub('PENDING')
        sh 'env | sort'
      }
    }
    
    stage('build-local') {
      when {
        not { expression { isBinaryBuild() } }
        expression { isUnix() }
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
      when { expression { isBinaryBuild() } }
      environment { PARENT_STAGE_NAME = "${STAGE_NAME}" }
      
      parallel {
        stage('notify') { steps { notifyGitHub('PENDING') } }
        stage('mini')   { steps { sh "bash ci_support/package.sh " + '${WORKSPACE} ${WORKSPACE}/ci_support/constructor-${STAGE_NAME}/' } }
        stage('huge')   { steps { sh "bash ci_support/package.sh " + '${WORKSPACE} ${WORKSPACE}/ci_support/constructor-${STAGE_NAME}/' } }
      }
    }
    
    stage('test-package') {
      when { expression { isBinaryBuild() } }
      environment { PARENT_STAGE_NAME = "${STAGE_NAME}" }
      
      parallel {
        stage('notify') { steps { notifyGitHub('PENDING') } }
        stage('mini')   { steps { testPackage(binary_size_suffix[STAGE_NAME], STAGE_NAME) } }
        stage('huge')   { steps { testPackage(binary_size_suffix[STAGE_NAME], STAGE_NAME) } }
      }
    }
    
    stage('deploy') {
      when { expression { isBinaryBuild() } }
      environment { PARENT_STAGE_NAME = "${STAGE_NAME}" }

      parallel {
        stage('notify') { steps { notifyGitHub('PENDING') } }
        stage('mini')   { steps { deployPackage(binary_size_suffix[STAGE_NAME]) } }
        stage('huge')   { steps { deployPackage(binary_size_suffix[STAGE_NAME]) } }
      }
    }
  }
  
  post {
    always {
      notifyGitHub("${currentBuild.result}")
      notifyEmail()
    }
  }
}
