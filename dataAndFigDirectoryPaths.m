function paths = dataAndFigDirectoryPaths(data_path)
    paths = struct();
    paths.data = [data_path, filesep];
    paths.individualRecordings = [paths.data, 'individualRecordings', filesep];
    paths.figures = [paths.data, 'figureExport', filesep];
    mkdir(paths.figures)
    paths.diffusionData = [paths.data, 'diffusion', filesep];
end