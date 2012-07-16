
; Copyright (C) 2012 Victor Semionov
; All rights reserved.
; 
; Redistribution and use in source and binary forms, with or without
; modification, are permitted provided that the following conditions are met:
;  * Redistributions of source code must retain the above copyright notice, this
;    list of conditions and the following disclaimer.
;  * Redistributions in binary form must reproduce the above copyright notice,
;    this list of conditions and the following disclaimer in the documentation
;    and/or other materials provided with the distribution.
; 
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
; ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
; DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
; FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
; SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
; CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
; OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#define MyAppName "NPAmp"
#define MyGUIAppName "GNPAmp"
; the version is also specified below for VersionInfoVersion
#define MyAppVersion "1.2.18"
#define MyAppPublisher "Victor Semionov"
#define MyPublisherURL "http://www.vsemionov.org/"
#define MyAppURL "http://www.vsemionov.org/npamp/"
#define MyAppExeName "npamp.exe"
#define MyGUIAppExeName "gnpamp.exe"
#define MyCopyrightHolder MyAppPublisher
#define MyCopyrightPeriod "2012"

[Setup]
AppId={{924BF823-C44E-4700-BE6C-D88EF72C00F6}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyPublisherURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={pf}\{#MyAppName}
DefaultGroupName={#MyAppName}
LicenseFile=..\LICENSE
OutputDir=.
OutputBaseFilename={#MyAppName} Setup {#MyAppVersion}
Compression=lzma
SolidCompression=yes

ChangesAssociations=yes

VersionInfoVersion=1.2.18.0
VersionInfoCopyright=Copyright (C) {#MyCopyrightPeriod} {#MyCopyrightHolder}

SetupIconFile=setup.ico
UninstallDisplayIcon={app}\{#MyGUIAppExeName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "fileassoc"; Description: "&Register file associations"; GroupDescription: "Other tasks:"

[Files]
Source: "..\vcredist_x86.exe"; DestDir: "{tmp}"; Flags: ignoreversion
Source: "..\dist\*.*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\README"; DestDir: "{app}"; DestName: "README.txt"; Flags: ignoreversion
Source: "..\AUTHORS"; DestDir: "{app}"; DestName: "AUTHORS.txt"; Flags: ignoreversion
Source: "..\LICENSE"; DestDir: "{app}"; DestName: "LICENSE.txt"; Flags: ignoreversion
Source: "..\CHANGES"; DestDir: "{app}"; DestName: "CHANGES.txt"; Flags: ignoreversion

[Icons]
Name: "{group}\{#MyGUIAppName}"; Filename: "{app}\{#MyGUIAppExeName}"
Name: "{commondesktop}\{#MyGUIAppName}"; Filename: "{app}\{#MyGUIAppExeName}"; Tasks: desktopicon
Name: "{group}\README"; Filename: "{app}\README.txt"

[Run]
Filename: "{tmp}\vcredist_x86.exe"; Parameters: "/qb!"
Filename: "{app}\{#MyGUIAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyGUIAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent
Filename: "{app}\README.txt"; Description: "View README"; Flags: shellexec nowait postinstall skipifdoesntexist skipifsilent unchecked

[Registry]
Root: HKCR; Subkey: ".npc"; ValueType: string; ValueName: ""; ValueData: "NPAmpFile"; Flags: uninsdeletevalue; Tasks: fileassoc
Root: HKCR; Subkey: "NPAmpFile"; ValueType: string; ValueName: ""; ValueData: "NPAmp File"; Flags: uninsdeletekey; Tasks: fileassoc
Root: HKCR; Subkey: "NPAmpFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\{#MyGUIAppExeName},1"; Tasks: fileassoc
Root: HKCR; Subkey: "NPAmpFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#MyGUIAppExeName}"" ""%1"""; Tasks: fileassoc
