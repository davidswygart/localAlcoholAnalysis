function paths = dataAndFigDirectoryPaths()
    mfilePath = [mfilename('fullpath'), '.m'];
    if contains(mfilePath,'LiveEditorEvaluationHelper')
        mfilePath = matlab.desktop.editor.getActiveFilename;
    end
    
    rootDir = dir(mfilePath).folder;

    paths = struct();
    paths.root = [rootDir, filesep];
    paths.data = [paths.root, 'data', filesep];
    paths.individualRecordings = [paths.data, 'individualRecordings', filesep];
    paths.figures = [paths.root, 'figureExport', filesep];
end