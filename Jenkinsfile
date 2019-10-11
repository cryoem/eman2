binary_size_suffix = ['mini':'', 'huge':'_huge']

def convertToNativePath(path) {
    return sh(returnStdout: true, script: "python -c 'import os; print(os.path.normpath(\"" + path + "\"))'").trim()
}

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
        switch(status) {
            case 'PENDING':
                message = 'Stage: ' + (env.PARENT_STAGE_NAME ?: STAGE_NAME)
                break
            case 'SUCCESS':
                message = 'Build succeeded!'
                break
            case 'FAILURE':
                message = 'Build failed!'
                break
            case 'ABORTED':
                message = 'Build aborted!'
                status == 'ERROR'
                break
        }

        context = "JenkinsCI/${AGENT_OS_NAME.capitalize()}"

        step([$class: 'GitHubCommitStatusSetter',
              contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: context],
              statusResultSource: [$class: 'ConditionalStatusResultSource',
                                   results: [[$class: 'AnyBuildResult', message: message, state: status]]]])
    }
}

def notifyEmail() {
    from    = "JenkinsCI ($AGENT_OS_NAME) <jenkins@jenkins>"
    body    = '''${SCRIPT, template="groovy-text.template"}'''
    subject = '$BUILD_STATUS! ' + "($GIT_BRANCH_SHORT - ${GIT_COMMIT_SHORT})" + ' #$BUILD_NUMBER'

    if(JOB_TYPE == "push" || NOTIFY_EMAIL == "true") {
        to      = "$GIT_AUTHOR_EMAIL"
    }
    if(JOB_TYPE == "cron") {
        to      = '$DEFAULT_RECIPIENTS'
        subject = '[cron] - ' + subject
    }

    emailext(to:        to,
             from:      from,
             subject:   subject,
             body:      body,
             attachLog: true
             )
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

// Can't be called to set pipeline envvars, because depends on CI_BUILD
def isContinuousBuild() {
    return (CI_BUILD == "1" && isMasterBranch()) || isReleaseBranch() || JOB_TYPE == "cron"
}

// Can't be called to set pipeline envvars, because depends on CI_BUILD
def isExperimentalBuild() {
    return CI_BUILD == "1" && !(isMasterBranch() || isReleaseBranch())
}

// Can't be called to set pipeline envvars, because depends on CI_BUILD indirectly
def isBinaryBuild() {
    return isContinuousBuild() || isExperimentalBuild()
}

def getBuildStabilityType() {
    if(isContinuousBuild())        return 'unstable'
    else if(isExperimentalBuild()) return 'experimental'
    else                           return 'NONE'
}

def getInstallerExt() {
    if(isUnix()) return 'sh'
    else         return 'exe'
}

def getDeployFileName(size_type) {
    stability_type = getBuildStabilityType()
    installer_ext  = getInstallerExt()

    return "eman2"+ binary_size_suffix[size_type] + ".${AGENT_OS_NAME}." + stability_type + "." + installer_ext
}

def testPackage(installer_file, installation_dir) {
    def script_file_base = convertToNativePath("tests/test_binary_installation.")
    installer_file       = convertToNativePath(installer_file)
    installation_dir     = convertToNativePath(installation_dir)
    
    if(isUnix())
        sh  "bash " + script_file_base + "sh "  + installer_file + " " + installation_dir
    else
        bat "call " + script_file_base + "bat " + installer_file + " " + installation_dir
}

def deployPackage(size_type='') {
    stability_type = getBuildStabilityType()
    installer_ext  = getInstallerExt()

    def sourceFile  = "eman2" + binary_size_suffix[size_type] + ".${AGENT_OS_NAME}." + installer_ext
    def targetFile  = getDeployFileName(size_type)
    def cdCommand   = "cd ${DEPLOY_PATH}/" + stability_type
    def mvCommand   = "mv " + sourceFile + " " + targetFile
    def execCommand = cdCommand + " && " + mvCommand
    
    sshPublisher(publishers: [
                              sshPublisherDesc(configName: 'Installer-Server',
                                               transfers:
                                                          [sshTransfer(sourceFiles:        sourceFile,
                                                                       removePrefix:       "",
                                                                       remoteDirectory:    stability_type,
                                                                       remoteDirectorySDF: false,
                                                                       cleanRemote:        false,
                                                                       excludes:           '',
                                                                       execCommand:        execCommand,
                                                                       execTimeout:        120000,
                                                                       flatten:            false,
                                                                       makeEmptyDirs:      false,
                                                                       noDefaultExcludes:  false,
                                                                       patternSeparator:   '[, ]+'
                                                                      )
                                                          ],
                                                          usePromotionTimestamp:   false,
                                                          useWorkspaceInPromotion: false,
                                                          verbose:                 true
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
    INSTALLERS_DIR = convertToNativePath("${HOME_DIR}/workspace/jenkins-eman-installers")

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
      environment {
        PARENT_STAGE_NAME = "${STAGE_NAME}"
        INSTALLER_EXT     = getInstallerExt()
      }
      
      parallel {
        stage('notify') { steps { notifyGitHub('PENDING') } }
        stage('mini')   { steps { testPackage("${WORKSPACE}/eman2" + binary_size_suffix[STAGE_NAME] + ".${AGENT_OS_NAME}.${INSTALLER_EXT}", "${INSTALLERS_DIR}/" + STAGE_NAME) } }
        stage('huge')   { steps { testPackage("${WORKSPACE}/eman2" + binary_size_suffix[STAGE_NAME] + ".${AGENT_OS_NAME}.${INSTALLER_EXT}", "${INSTALLERS_DIR}/" + STAGE_NAME) } }
      }
    }
    
    stage('deploy') {
      when { expression { isBinaryBuild() } }
      environment { PARENT_STAGE_NAME = "${STAGE_NAME}" }

      parallel {
        stage('notify') { steps { notifyGitHub('PENDING') } }
        stage('mini')   { steps { deployPackage(STAGE_NAME) } }
        stage('huge')   { steps { deployPackage(STAGE_NAME) } }
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
