#Requires -Version 7
param([switch] $Start, [switch] $Install)

$ErrorActionPreference = 'Stop'
Push-Location (Split-Path -Parent $PSScriptRoot)
try {
    if (!(Test-Path -LiteralPath '.venv/Scripts/python.exe')) {
        python -m venv .venv
        if ($LASTEXITCODE -ne 0) { throw '创建本地虚拟环境失败' }
        $Install = $true
    }
    if ($Start -and !$Install) {
        & '.venv/Scripts/python.exe' -c 'import streamlit, plotly, aphcp_control' 2>$null
        if ($LASTEXITCODE -ne 0) { $Install = $true }
    }
    if (!$Start -or $Install) {
        & '.venv/Scripts/python.exe' -m pip install -r requirements-local.txt
        if ($LASTEXITCODE -ne 0) { throw '安装管理端依赖失败' }
        & '.venv/Scripts/python.exe' -m pip install --no-deps -e .
        if ($LASTEXITCODE -ne 0) { throw '安装管理端包失败' }
    }
    if ($Start) {
        & '.venv/Scripts/python.exe' -m streamlit run src/aphcp_dashboard/app.py --server.address 127.0.0.1 --server.port 8501
        if ($LASTEXITCODE -ne 0) { throw '启动 Dashboard 失败' }
    }
} finally {
    Pop-Location
}
