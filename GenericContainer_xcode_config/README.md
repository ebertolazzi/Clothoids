## xcconfig

A pretty good project starts with a pretty good configuration. Really, the more constraints you enable, the more stable and higher quality your project will be for it.

These `.xcconfig` files define build configurations for Xcode without having to go into the Xcode interface itself to make changes. Also, they can be easily shared across projects.

### Instructions

Projects may be built with separate schemes: Debug, Release, and Test.

Each scheme (Debug/Release/Test) specifies it's own set of criteria. Criteria defined in `Base.xcconfig` is shared and inherited by the other files.

Check out a few of my other projects for an [example](https://github.com/piemonte/PBJVision/).

### Usage

As of Xcode 7.0, these `xcconfig` files can be setup for a project by going to the Project Editor's 'Info' tab, then 'Configurations'.

If you would like to adopt these for an existing project, I recommend commenting out each line and slowly enabling them as you fix and update the code.

### References

[Project Editor Help: Add a Build Configuration](https://developer.apple.com/library/ios/recipes/xcode_help-project_editor/Articles/BasingBuildConfigurationsonConfigurationFiles.html)
