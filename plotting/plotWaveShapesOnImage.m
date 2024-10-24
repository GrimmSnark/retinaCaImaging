function plotWaveShapesOnImage(imageStack, waveTable, frame)

imshow(imageStack(:,:,frame));
hold on

waveTableIndx = waveTable.Frame == frame;
object2Plot = waveTable(waveTableIndx,:);

for i = 1:height(object2Plot)
rectangle('Position', object2Plot.BoundingBox(i,:), 'EdgeColor', 'r');
scatter(object2Plot.Centroid(i,1), object2Plot.Centroid(i,2));
end

end